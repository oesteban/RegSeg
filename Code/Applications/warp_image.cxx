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
	std::vector<size_t> grid_size;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("images,I", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("surfaces,S", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken(),	"moving image file")
			("displacement-field,F", bpo::value < std::string >(&fieldname), "forward displacement field" )
			("inverse-displacement-field,R", bpo::value < std::string >(&invfieldname), "backward displacement field" )
			("mask,M", bpo::value< std::string >(&maskfile), "mask file" )
			//("compute-inverse", bpo::bool_switch(), "compute precise inversion of the input field (requires -F)")
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("mask-inputs", bpo::bool_switch(), "use deformed mask to filter input files")
			("grid-size,g", bpo::value< std::vector<size_t> >(&grid_size)->multitoken(), "size of grid of bspline control points (default is 10x10x10)");

	bpo::variables_map vm;
	bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
	bpo::notify(vm);

	if (vm.count("help") || vm.size() == 0 ) {
		std::cout << all_desc << std::endl;
		return 1;
	}

	if ( !vm.count("displacement-field") && !vm.count("inverse-displacement-field") ) {
		std::cerr << "One of -F or -R options should be specified" << std::endl;
		return 1;
	}

	if ( vm.count("displacement-field") && vm.count("inverse-displacement-field") ) {
		std::cerr << "-F or -R options are mutually exclusive" << std::endl;
		return 1;
	}

	bool isFwdField = vm.count("displacement-field");
	bool computeInv = vm.count("compute-inverse");

	DisplacementFieldPointer field, field_inv, input_field;
	DisplacementFieldReaderPointer fread = DisplacementFieldReaderType::New();
	fread->SetFileName( isFwdField?fieldname:invfieldname );
	fread->Update();
	input_field = fread->GetOutput();

	const VectorType* ofb;
	VectorType* ifb;

	if( isFwdField ) {
		field = input_field;
		field_inv = FieldType::New();
		field_inv->SetRegions( field->GetLargestPossibleRegion());
		field_inv->SetOrigin( field->GetOrigin());
		field_inv->SetDirection( field->GetDirection() );
		field_inv->Allocate();

		ofb = field->GetBufferPointer();
		ifb = field_inv->GetBufferPointer();
	} else {
		field_inv = input_field;
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

		if (mask.IsNotNull() && (vm.count("mask-inputs") && vm["mask-inputs"].as<bool>() ) ) {
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
	typename FieldType::SizeType size;
	if( vm.count("grid-size") ){
		if( grid_size.size() == 1 ) {
			size.Fill( grid_size[0] );
		}
		else if ( grid_size.size() == FieldType::ImageDimension ) {
			for( size_t i = 0; i < FieldType::ImageDimension; i++)
				size[i] = grid_size[i];
		}
		else {
			std::cout << "error with grid size" << std::endl;
			return 1;
		}
	}

	TPointer tf_mesh = Transform::New();
	tf_mesh->SetControlPointsSize(size);
	tf_mesh->SetPhysicalDomainInformation( field );
	tf_mesh->SetField( field );
	tf_mesh->ComputeCoefficients();

	// TransformPointer tf_inv = TransformType::New();
	// tf_inv->SetDisplacementField( field_inv );


	for( size_t i = 0; i<movingSurfaceNames.size(); i++){
		MeshReaderPointer r = MeshReaderType::New();
		r->SetFileName( movingSurfaceNames[i] );
		r->Update();

		MeshPointer cur_mesh = r->GetOutput();
		PointsIterator p_it = cur_mesh->GetPoints()->Begin();
		PointsIterator p_end = cur_mesh->GetPoints()->End();

		MeshPointType p;
		while ( p_it!=p_end ) {
			tf_mesh->AddOffGridPos(p_it.Value());
			++p_it;
		}
	}

	tf_mesh->Interpolate();

	size_t pointId = 0;
	for( size_t i = 0; i<movingSurfaceNames.size(); i++){
		MeshReaderPointer r = MeshReaderType::New();
		r->SetFileName( movingSurfaceNames[i] );
		r->Update();

		MeshPointer cur_mesh = r->GetOutput();
		PointsIterator p_it = cur_mesh->GetPoints()->Begin();
		PointsIterator p_end = cur_mesh->GetPoints()->End();

		MeshPointType p;
		while ( p_it!=p_end ) {
			p = p_it.Value();
			p_it.Value() += tf_mesh->GetOffGridValue(pointId);
			++p_it;
		}

		MeshWriterPointer wmesh = MeshWriterType::New();
		std::stringstream ss;
		ss << outPrefix << "_warped_" << i << ".vtk";
		wmesh->SetFileName( ss.str().c_str() );
		wmesh->SetInput( cur_mesh );
		wmesh->Update();
	}
}

