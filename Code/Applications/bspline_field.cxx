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

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <sstream>


int main(int argc, char *argv[]) {
	std::string outPrefix = "";
	std::string maskfile;
	std::vector< std::string > fixedImageNames, movingSurfaceNames,coefficientImageNames;
	std::vector<size_t> grid_size;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("coeff-images,C", bpo::value < std::vector<std::string> > (&coefficientImageNames )->multitoken(), "coefficient image(s)" )
			("images,I", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("surfaces,S", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken(),	"moving image file")
			("mask,M", bpo::value< std::string >(&maskfile), "mask file" )
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("num-threads", bpo::value < unsigned int >()->default_value(NUM_THREADS), "use num-threads")
			("grid-size,g", bpo::value< std::vector<size_t> >(&grid_size)->multitoken(), "size of grid of bspline control points (default is 10x10x10)")
			("mask-inputs", bpo::bool_switch(), "use deformed mask to filter input files");

	bpo::variables_map vm;

	try {
		bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);

		if (vm.count("help")) {
			std::cout << all_desc << std::endl;
			return 1;
		}

		bpo::notify(vm);
	} catch ( bpo::error& e ) {
		std::cerr << "Error: " << e.what() << std::endl << std::endl;
		std::cerr << all_desc << std::endl;
		return 1;
	}

	ReaderPointer readref = ReaderType::New();
	readref->SetFileName( fixedImageNames[0] );
	readref->Update();
	ChannelPointer ref = readref->GetOutput();
	typename ChannelType::DirectionType ref_dir(ref->GetDirection());
	typename ChannelType::PointType ref_orig = ref->GetOrigin();

	typename ChannelType::DirectionType itk;
	itk.SetIdentity();
	itk(0,0)=-1.0;
	itk(1,1)=-1.0;

	typename ChannelType::DirectionType int_dir(itk * ref_dir);
	typename ChannelType::PointType int_orig( itk * ref_orig );
	ref->SetDirection( int_dir );
	ref->SetOrigin( int_orig );

	typename CoefficientsType::SizeType size;
	typename CoefficientsType::SpacingType spacing;

	size.Fill(10);
	CoefficientsImageArray coeffs;

	TPointer transform = Transform::New();

#ifndef NDEBUG
	transform->SetNumberOfThreads( 2 );
#endif

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

		transform->SetControlGridSize( size );
		transform->SetDomainExtent(ref);
		transform->SetOutputReference( ref );
		transform->Initialize();

		coeffs = transform->GetCoefficientsImages();
		spacing = coeffs[0]->GetSpacing();
		size_t numPix = coeffs[0]->GetLargestPossibleRegion().GetNumberOfPixels();

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
#ifdef NDEBUG
			rndit.ReinitializeSeed();
#endif
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

			// Remove median
			SubtractPointer sub = SubtractFilter::New();
			sub->SetInput1( coeffs[i] );
			sub->SetConstant( sample[ static_cast<size_t>( floor( 0.5 * numPix + 0.5f ) ) ] );
			sub->Update();
			coeffs[i] = sub->GetOutput();

			// Compute maximum displacement (check both endpoints)
			MaxCalcPointer max = MaxCalc::New();
			max->SetImage( coeffs[i] );
			max->Compute();

			float immax = max->GetMaximum();
			if ( fabs( max->GetMinimum() ) > fabs(immax) ) {
				immax = fabs( max->GetMinimum() );
			}

			// Scale to maximum diffeomorphic displacement
			float max_disp = (spacing[i] * 0.39);

			MultiplyPointer m = MultiplyFilter::New();
			m->SetInput1( coeffs[i] );
			m->SetConstant( max_disp / immax );
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

	typename FieldType::Pointer field = transform->GetDisplacementField();
	typename FieldWriter::Pointer ff = FieldWriter::New();
	ff->SetInput( field );
	ff->SetFileName( (outPrefix + "_field.nii.gz").c_str() );
	ff->Update();

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
	DisplacementFieldTransformPointer tf_inv = DisplacementFieldTransformType::New();
	tf_inv->SetDisplacementField( transform->GetInverseDisplacementField() );

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
		ss.str("");
		ss << outPrefix << "_warped_" << i << ".vtk";
		wmesh->SetFileName( ss.str().c_str() );
		wmesh->SetInput( mesh );
		wmesh->Update();
	}
}
