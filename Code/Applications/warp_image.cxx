/*
 * warp_image.cxx
 *
 *  Created on: Mar 27, 2014
 *      Author: oesteban
 */

#include "warp_image.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <sstream>

int main(int argc, char *argv[]) {
	std::string outPrefix = "displ";
	std::string fieldname,maskfile,invfieldname;
	std::vector< std::string > fixedImageNames, movingSurfaceNames;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("images,I", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("surfaces,S", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken(),	"moving image file")
			("displacement-field,F", bpo::value < std::string >(&fieldname), "forward displacement field" )
			("inverse-displacement-field,R", bpo::value < std::string >(&invfieldname), "backward displacement field" )
			("mask,M", bpo::value< std::string >(&maskfile), "mask file" )
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files");

	bpo::variables_map vm;
	bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
	bpo::notify(vm);

	if (vm.count("help") || vm.size() == 0 ) {
		std::cout << all_desc << std::endl;
		return 1;
	}

	if (vm.count("displacement-field")==0 && vm.count("inverse-displacement-field")==0 ) {
		std::cerr << "One of -F or -R options should be specified" << std::endl;
		return 1;
	}

	DisplacementFieldReaderPointer fread = DisplacementFieldReaderType::New();
	DisplacementFieldPointer field, field_inv;

	const VectorType* ofb;
	VectorType* ifb;

	if( vm.count("displacement-field") ) {
		fread->SetFileName( fieldname );
		fread->Update();
		field = fread->GetOutput();
		field_inv = FieldType::New();
		field_inv->SetRegions( field->GetLargestPossibleRegion());
		field_inv->SetOrigin( field->GetOrigin());
		field_inv->SetDirection( field->GetDirection() );
		field_inv->Allocate();

		ofb = field->GetBufferPointer();
		ifb = field_inv->GetBufferPointer();
	} else {
		fread->SetFileName( invfieldname );
		fread->Update();
		field_inv = fread->GetOutput();
		field = FieldType::New();
		field->SetRegions( field_inv->GetLargestPossibleRegion());
		field->SetOrigin( field_inv->GetOrigin());
		field->SetDirection( field_inv->GetDirection() );
		field->Allocate();

		ofb = field_inv->GetBufferPointer();
		ifb = field->GetBufferPointer();
	}

	size_t nPix = field->GetLargestPossibleRegion().GetNumberOfPixels();
	VectorType v;
	for( size_t i = 0; i<nPix; i++ ) {
		v = *( ofb+i );
		*( ifb + i ) = -v;
	}

	// Direction issues
	typename ChannelType::DirectionType itk;
	itk.SetIdentity();
	itk(0,0)=-1.0;
	itk(1,1)=-1.0;

	typename ReaderType::Pointer rref = ReaderType::New();
	rref->SetFileName( fixedImageNames[0] );
	rref->Update();
	typename ChannelType::Pointer ref = rref->GetOutput();

	typename ChannelType::DirectionType ref_dir = ref->GetDirection();
	typename ChannelType::PointType ref_orig = ref->GetOrigin();
	typename ChannelType::DirectionType int_dir(itk * ref_dir);
	typename ChannelType::PointType int_orig( itk * ref_orig );

	TransformPointer transform = TransformType::New();
	transform->SetDisplacementField(field);

	// Read and transform mask, if present
	MaskPointer mask;
	if (vm.count( "mask" ) ) {
		typename ReaderType::Pointer rmask = ReaderType::New();
		rmask->SetFileName( maskfile );
		rmask->Update();

		typename ChannelType::Pointer im = rmask->GetOutput();
		im->SetDirection( int_dir );
		im->SetOrigin( int_orig );

		typename Binarize::Pointer bin = Binarize::New();
		bin->SetInput( im );
		bin->SetLowerThreshold( 0.01 );
		bin->SetOutsideValue( 0 );
		bin->SetInsideValue( 1 );
		bin->Update();

		WarpMaskFilterPointer wrp = WarpMaskFilter::New();
		wrp->SetInterpolator( WarpMaskInterpolator::New() );
		wrp->SetOutputParametersFromImage( bin->GetOutput() );
		wrp->SetInput( bin->GetOutput() );
		wrp->SetDisplacementField( field );
		wrp->Update();
		mask = wrp->GetOutput();

		mask->SetDirection( ref_dir );
		mask->SetOrigin( ref_orig );

		typename MaskWriter::Pointer wm = MaskWriter::New();
		wm->SetInput( mask );
		wm->SetFileName( (outPrefix + "_mask_warped.nii.gz").c_str() );
		wm->Update();
	}

	// Read and transform images if present
	for( size_t i = 0; i<fixedImageNames.size(); i++) {
		std::stringstream ss;
		typename WriterType::Pointer w = WriterType::New();

		typename ReaderType::Pointer r = ReaderType::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();

		typename ChannelType::Pointer im = r->GetOutput();
		im->SetDirection( int_dir );
		im->SetOrigin( int_orig );

		WarpFilterPointer wrp = WarpFilter::New();
		wrp->SetInterpolator( WarpInterpolator::New() );
		wrp->SetOutputParametersFromImage( im );
		wrp->SetInput( im );
		wrp->SetDisplacementField( field );
		wrp->Update();

		ThresholdPointer th = ThresholdFilter::New();
		th->SetInput( wrp->GetOutput() );
		th->ThresholdBelow( 0.0 );
		th->SetOutsideValue( 0.0 );
		th->Update();

		typename ChannelType::Pointer im_wrp = th->GetOutput();
		im_wrp->SetDirection( ref_dir );
		im_wrp->SetOrigin( ref_orig );

		if (mask.IsNotNull()) {
			typename MaskFilter::Pointer mm = MaskFilter::New();
			mm->SetMaskImage( mask );
			mm->SetInput( im_wrp );
			mm->Update();
			im_wrp = mm->GetOutput();
		}

		ss.str("");
		ss << outPrefix << "_warped_" << i << ".nii.gz";
		w->SetInput( im_wrp );
		w->SetFileName( ss.str().c_str() );
		w->Update();
	}


	// Warp surfaces --------------------------------------------------
	TransformPointer tf_inv = TransformType::New();
	tf_inv->SetDisplacementField( field_inv );


	for( size_t i = 0; i<movingSurfaceNames.size(); i++){
		MeshReaderPointer r = MeshReaderType::New();
		r->SetFileName( movingSurfaceNames[i] );
		r->Update();

		MeshPointer mesh = r->GetOutput();

		PointsIterator p_it = mesh->GetPoints()->Begin();
		PointsIterator p_end = mesh->GetPoints()->End();

		MeshPointType p;
		while ( p_it!=p_end ) {
			p = p_it.Value();
			p_it.Value() = tf_inv->TransformPoint( p );
			++p_it;
		}

		MeshWriterPointer wmesh = MeshWriterType::New();
		std::stringstream ss;
		ss << outPrefix << "_warped_" << i << ".vtk";
		wmesh->SetFileName( ss.str().c_str() );
		wmesh->SetInput( mesh );
		wmesh->Update();
	}
}

