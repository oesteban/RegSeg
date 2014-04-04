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
	std::string fieldname,maskfile;
	std::vector< std::string > fixedImageNames, movingSurfaceNames;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("images,I", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("surfaces,S", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken(),	"moving image file")
			("displacement-field,F", bpo::value < std::string >(&fieldname), "displacement field" )
			("mask,M", bpo::value< std::string >(&maskfile), "mask file" )
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files");

	bpo::variables_map vm;
	bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
	bpo::notify(vm);

	if (vm.count("help") || vm.size() == 0 ) {
		std::cout << all_desc << std::endl;
		return 1;
	}

	DisplacementFieldReaderPointer fread = DisplacementFieldReaderType::New();
	fread->SetFileName( fieldname );
	fread->Update();
	DisplacementFieldPointer field = fread->GetOutput();

	// Direction issues
	typename ChannelType::DirectionType itk;
	itk.SetIdentity();
	itk(0,0)=-1.0;
	itk(1,1)=-1.0;

	typename ReaderType::Pointer rref = ReaderType::New();
	rref->SetFileName( fixedImageNames[0] );
	rref->Update();

	typename ChannelType::Pointer im = rref->GetOutput();
	typename ChannelType::DirectionType dir = im->GetDirection();
	typename ChannelType::PointType ref_orig = im->GetOrigin();


	// Read and transform mask, if present
	MaskPointer mask;

	if (vm.count( "mask" ) ) {
		typename ReaderType::Pointer rmask = ReaderType::New();
		rmask->SetFileName( maskfile );
		rmask->Update();

		ChannelPointer im = rmask->GetOutput();
		im->SetDirection( dir * itk );
		im->SetOrigin( itk * ref_orig );

		WarpFilterPointer res = WarpFilter::New();
		res->SetInput( im );
		res->SetOutputParametersFromImage( rmask->GetOutput() );
		res->SetInterpolator( NearestNeighborInterpolateImageFunction::New() );
		res->SetDisplacementField( field );
		res->Update();

		ChannelPointer im_res = res->GetOutput();
		im_res->SetDirection( dir );
		im_res->SetOrigin( ref_orig );

		typename Binarize::Pointer bin = Binarize::New();
		bin->SetInput( im_res );
		bin->SetLowerThreshold( 0.01 );
		bin->SetOutsideValue( 0 );
		bin->SetInsideValue( 1 );
		bin->Update();
		mask = bin->GetOutput();

		typename MaskWriter::Pointer wm = MaskWriter::New();
		wm->SetInput( mask );
		wm->SetFileName( (outPrefix + "_mask_warped.nii.gz").c_str() );
		wm->Update();
	}

	// Read and transform images if present
	for( size_t i = 0; i<fixedImageNames.size(); i++) {
		typename ReaderType::Pointer r = ReaderType::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();

		typename ChannelType::Pointer im = r->GetOutput();
		typename ChannelType::DirectionType dir = im->GetDirection();
		typename ChannelType::PointType ref_orig = im->GetOrigin();

		im->SetDirection( dir * itk );
		im->SetOrigin( itk * ref_orig );

		WarpFilterPointer res = WarpFilter::New();
		res->SetInterpolator( itk::BSplineInterpolateImageFunction< ChannelType, ScalarType >::New() );
		res->SetOutputParametersFromImage( im );
		res->SetInput( im );
		res->SetDisplacementField( field );
		res->Update();

		typename ChannelType::Pointer im_res = res->GetOutput();
		im_res->SetDirection( dir );
		im_res->SetOrigin( ref_orig );

		if (mask.IsNotNull()) {
			typename MaskFilter::Pointer mm = MaskFilter::New();
			mm->SetMaskImage( mask );
			mm->SetInput( im_res );
			mm->Update();

			im_res = mm->GetOutput();
		}

		std::stringstream ss;
		ss.str("");
		ss << outPrefix << "_warped_" << i << ".nii.gz";
		typename WriterType::Pointer w = WriterType::New();
		w->SetInput( im_res );
		w->SetFileName( ss.str().c_str() );
		w->Update();

	}

	TransformPointer transform;
	if ( movingSurfaceNames.size() > 0 ) {
		transform = TransformType::New();
		transform->SetDisplacementField(field);
	}

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
			p_it.Value()+= p - transform->TransformPoint( p );
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

