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

#include <sstream>

int main(int argc, char *argv[]) {
	std::string outPrefix = "";
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
		transform->SetPhysicalDomainInformation( ref );
		transform->SetOutputReference( ref );
		transform->UpdateField();
		coeffs = transform->GetCoefficientsImages();
		spacing = coeffs[0]->GetSpacing();
		size_t numPix = coeffs[0]->GetLargestPossibleRegion().GetNumberOfPixels();

		std::cout << size << std::endl;
		std::cout << spacing << std::endl;

		typename CoefficientsType::RegionType region;
		typename CoefficientsType::IndexType start;
		typename CoefficientsType::SizeType regionSize;

		for( size_t i = 0; i<DIMENSION; i++ ) {
			start[i] = static_cast<size_t>( floor( size[i] * 0.35 + 0.5f ) ) -1;
			regionSize[i] = size[i] - 2.0 * start[i];
		}
		region.SetIndex( start );
		region.SetSize( regionSize );

		for( size_t i = 0; i< DIMENSION; i++) {
			RandomIterator rndit( coeffs[i], region );
			rndit.SetNumberOfSamples( static_cast<size_t>( floor( numPix * 0.08 + 0.5f ) ) );

			for(rndit.GoToBegin(); !rndit.IsAtEnd(); ++rndit){
				float r = -1.0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/2.0));
				rndit.Set( r );
			}

			SigmaArrayType sigma;
			sigma.Fill(12.0);

			SmoothingFilterPointer s = SmoothingFilterType::New();
			s->SetInput( coeffs[i] );
			s->SetSigmaArray( sigma );
			s->Update();
			coeffs[i] = s->GetOutput();

			ScalarType *buff = coeffs[i]->GetBufferPointer();
			std::vector< ScalarType > sample;

			for( size_t j = 0; j < numPix; j++ ) {
				sample.push_back( *(buff + j) );
			}

			std::sort( sample.begin(), sample.end() );

			SubtractPointer sub = SubtractFilter::New();
			sub->SetInput1( coeffs[i] );
			sub->SetConstant( sample[ static_cast<size_t>( floor( 0.5 * numPix + 0.5f ) ) ] );
			sub->Update();
			coeffs[i] = sub->GetOutput();

			MaxCalcPointer max = MaxCalc::New();
			max->SetImage( coeffs[i] );
			max->Compute();

			float immax = max->GetMaximum();
			if ( fabs( max->GetMinimum() ) > fabs(immax) ) {
				immax = fabs( max->GetMinimum() );
			}

			float scale = (spacing[i] * 0.35) / immax;

			MultiplyPointer m = MultiplyFilter::New();
			m->SetInput1( coeffs[i] );
			m->SetConstant( scale );
			m->Update();
			coeffs[i] = m->GetOutput();

			CoefficientsWriterPointer w = CoefficientsWriterType::New();
			std::stringstream ss;
			ss << outPrefix << "_coeffs_" << i << ".nii.gz";
			w->SetFileName( ss.str().c_str() );
			w->SetInput( coeffs[i] );
			w->Update();
		}
	}

	// Set coefficients
	transform->SetCoefficientsImages( coeffs );
	transform->UpdateField();

	typename ComponentsWriter::Pointer f = ComponentsWriter::New();
	std::stringstream ss;
	ss << outPrefix << "_field";
	f->SetFileName( ss.str().c_str() );
	f->SetInput( transform->GetField() );
	f->Update();


	transform->Interpolate();
	f->SetInput( transform->GetOutputField() );
	ss.str("");
	ss << outPrefix << "_field_hires";
	f->SetFileName( ss.str().c_str() );
	f->Update();


	// Read and transform images if present
	for( size_t i = 0; i<fixedImageNames.size(); i++) {
		ReaderPointer r = ReaderType::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();


		ResamplePointer res = ResampleFilter::New();
		res->SetInput( r->GetOutput() );
		res->SetReferenceImage( r->GetOutput() );
		res->SetUseReferenceImage(true);
		res->SetInterpolator( BSplineInterpolateImageFunction::New() );
		res->SetTransform( transform );
		res->Update();

		ss.str("");
		ss << outPrefix << "_resampled_" << i << ".nii.gz";
		WriterPointer w = WriterType::New();
		w->SetInput( res->GetOutput() );
		w->SetFileName( ss.str().c_str() );
		w->Update();
	}

}
