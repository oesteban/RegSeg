// --------------------------------------------------------------------------------------
// File:          bspline_field.cxx
// Date:          Mar 5, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "bspline_field.h"


int main(int argc, char *argv[]) {
	std::string outPrefix;
	std::vector< std::string > fixedImageNames, movingSurfaceNames,coefficientImageNames;
	std::vector<size_t> grid_size;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("coeff-images,C", bpo::value < std::vector<std::string> > (&coefficientImageNames )->multitoken(), "coefficient image(s)" )
			("target-images,T", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("moving-surfaces,M", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken(),	"moving image file")
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("grid-size,g", bpo::value< std::vector<size_t> >(&grid_size)->multitoken(), "size of grid of bspline control points (default is 10x10x10)");

	bpo::variables_map vm;
	bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
	bpo::notify(vm);

	if (vm.count("help") || vm.size() == 0 ) {
		std::cout << all_desc << std::endl;
		return 1;
	}

	ReaderPointer readref = ReaderType::New();
	readref->SetFileName( fixedImageNames[0] );
	readref->Update();
	ChannelPointer ref = readref->GetOutput();

	typename CoefficientsType::SizeType size;
	typename CoefficientsType::SpacingType spacing;
	size.Fill(10);
	CoefficientsImageArray coeffs;

	TPointer transform = Transform::New();
	transform->SetPhysicalDomainInformation( ref );
	transform->SetOutputReference( ref );

	if (vm.count("coeff-images")) {
		std::cout << "coefficient images mode not implemented" << std::endl;
		// set size
		return 0;
	} else {
		if( vm.count("grid-size") ){
			if( grid_size.size() == 1 ) {
				size.Fill( grid_size[0] );
			}
			else if ( grid_size.size() == DIMENSION ) {
				for( size_t i = 0; i < DIMENSION; i++) size[i] = grid_size[i];
			}
			else {
				std::cout << "error with grid size" << std::endl;
				return 1;
			}
		}

		transform->SetControlPointsSize( size );
		transform->UpdateField();
		coeffs = transform->GetCoefficientsImages();
		spacing = coeffs[0]->GetSpacing();

		std::cout << size << std::endl;
		std::cout << spacing << std::endl;

		for( size_t i = 0; i< DIMENSION; i++) {
			float scale = spacing[i] * 0.35;

			RandomSourcePointer src = RandomSourceType::New();
			src->SetMax( scale );
			src->SetMin( -1.0 * scale );
			src->SetSize( coeffs[i]->GetLargestPossibleRegion().GetSize() );
			src->SetSpacing( coeffs[i]->GetSpacing() );
			src->SetDirection( coeffs[i]->GetDirection() );
			src->SetOrigin( coeffs[i]->GetOrigin() );
			//src->SetScale( scale );
			src->Update();
			coeffs[i] = src->GetOutput();

			std::cout << scale << "mm." << std::endl;
			CoefficientsWriterPointer w = CoefficientsWriterType::New();
			w->SetFileName( "test.nii.gz" );
			w->SetInput( src->GetOutput() );
			w->Update();
		}
	}

	// Set coefficients




	// Read and transform images if present
	std::vector< ChannelPointer > images;
	for( size_t i = 0; i<fixedImageNames.size(); i++) {
		ReaderPointer r = ReaderType::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();
		images.push_back( r->GetOutput() );

	}

}
