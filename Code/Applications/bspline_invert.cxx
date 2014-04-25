// --------------------------------------------------------------------------------------
// File:          bspline_invert.cxx
// Date:          Apr 21, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//


#include "bspline_invert.h"

#include <vector>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <sstream>

int main(int argc, char *argv[]) {
	std::string outPrefix = "bsinv";
	std::string fixedname, fieldname;
	std::vector< std::string > coeffnames;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("coefficients-images,C", bpo::value < std::vector< std::string > >(&coeffnames)->multitoken()->required(), "forward bspline coefficients" )
			("displacement-field,D", bpo::value< std::string >(&fieldname)->required(), "displacement field" )
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files");

	bpo::variables_map vm;

	try {
		bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
		if (vm.count("help") || vm.size() == 0 ) {
			std::cout << all_desc << std::endl;
			return 1;
		}
		bpo::notify(vm);
	} catch ( bpo::error& e ) {
		std::cerr << "Error: " << e.what() << std::endl << std::endl;
		std::cerr << all_desc << std::endl;
		return 1;
	}

	DisplacementFieldReaderPointer rfield = DisplacementFieldReaderType::New();
	rfield->SetFileName( fieldname );
	rfield->Update();
	DisplacementFieldPointer field = rfield->GetOutput();

	CoefficientsImageArray coeffs;

	if ( coeffnames.size() == 3 ) {
		for( size_t i=0; i<coeffnames.size(); i++) {
			CoeffReaderPointer rcoeff = CoeffReaderType::New();
			rcoeff->SetFileName( coeffnames[i] );
			rcoeff->Update();
			coeffs[i] = rcoeff->GetOutput();
		}

	} else if ( coeffnames.size() == 1 ) {
		DisplacementFieldReaderPointer fread = DisplacementFieldReaderType::New();
		fread->SetFileName( coeffnames[0] );
		fread->Update();

		// TODO: split coefficients components
	}

	// Define a grid size
	typename CoefficientsType::SizeType size = coeffs[0]->GetLargestPossibleRegion().GetSize();
	typename CoefficientsType::SpacingType spacing = coeffs[0]->GetSpacing();

	// Set up a sparse matrix transform
	TPointer tf = Transform::New();
#ifndef NDEBUG
	tf->SetNumberOfThreads( 2 );
#endif
	tf->SetControlPointsSize( size );
	tf->SetPhysicalDomainInformation( field );

	VectorType disp;
	typename FieldType::IndexType idx;
	typename FieldType::PointType p,p_new;

	size_t nPix = field->GetLargestPossibleRegion().GetNumberOfPixels();
	tf->SetNumberOfSamples( nPix );
	const VectorType* fbuff = field->GetBufferPointer();

	for ( size_t i=0; i<nPix; i++) {
		// Find all new targets
		idx = field->ComputeIndex( i );
		field->TransformIndexToPhysicalPoint( idx, p );

		disp = *( fbuff + i );

		// Set inverse vector in new targets
		tf->SetOffGridPos( i, p + disp );
		tf->SetOffGridValue( i, -disp );
	}

	// Approximate coefficients
	tf->InvertField();


	typename FieldType::ConstPointer field_inv = tf->GetField();

	typename FieldWriter::Pointer ff = FieldWriter::New();
	ff->SetInput( field_inv );
	ff->SetFileName( (outPrefix + "_invfield.nii.gz").c_str() );
	ff->Update();
}
