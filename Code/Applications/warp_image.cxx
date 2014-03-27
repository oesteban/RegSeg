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
	std::string fieldname;
	std::vector< std::string > fixedImageNames;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("displacement-field,F", bpo::value < std::string >(&fieldname), "displacement field" )
			("images,I", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
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

	// Read and transform images if present
	for( size_t i = 0; i<fixedImageNames.size(); i++) {
		typename ReaderType::Pointer r = ReaderType::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();

		typename ChannelType::Pointer im = r->GetOutput();
		typename ChannelType::DirectionType dir = im->GetDirection();
		typename ChannelType::PointType ref_orig = im->GetOrigin();

		typename ChannelType::DirectionType itk;
		itk.SetIdentity();
		itk(0,0)=-1.0;
		itk(1,1)=-1.0;
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

		std::stringstream ss;
		ss.str("");
		ss << outPrefix << "_warped_" << i << ".nii.gz";
		typename WriterType::Pointer w = WriterType::New();
		w->SetInput( im_res );
		w->SetFileName( ss.str().c_str() );
		w->Update();

	}
}

