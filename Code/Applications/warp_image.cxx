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
#include <vector>

void conflicting_options(const boost::program_options::variables_map & vm,
						 const std::vector<std::string> & opts) {

	size_t optsset = 0;
	std::stringstream sopts;
	for (size_t i = 0; i<opts.size(); i++) {
		optsset+= int((vm.count(opts[i]) && !vm[opts[i]].defaulted()));
		sopts << " --" << opts[i];
		if (i!=(opts.size()-1))
			sopts << ",";
	}

	if (optsset != 1) {
		std::string s("One option of ");
		throw std::logic_error(std::string("One option of" + sopts.str() + " must be specified."));
	}
}

int main(int argc, char *argv[]) {
	std::string outPrefix = "displ";
	std::string maskfile;
	std::vector< std::string > fixedImageNames, movingSurfaceNames, fieldname, invfieldname;
	std::vector<size_t> grid_size;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("images,i", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("surfaces,s", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken(),	"moving image file")
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("mask,m", bpo::value< std::string >(&maskfile), "mask file" )
			("mask-inputs", bpo::bool_switch(), "use deformed mask to filter input files")
			("field,F", bpo::value < std::vector< std::string > >(&fieldname), "forward displacement field" )
			("inv-field,R", bpo::value < std::vector< std::string > >(&invfieldname), "backward displacement field" )
			("coeff,C", bpo::value < std::vector< std::string > >(&fieldname), "forward displacement field" )
			("inv-coeff,I", bpo::value < std::vector< std::string > >(&invfieldname), "backward displacement field" )
			//("compute-inverse", bpo::bool_switch(), "compute precise inversion of the input field (requires -F)")
			("grid-size,g", bpo::value< std::vector<size_t> >(&grid_size)->multitoken(), "size of grid of bspline control points (default is 10x10x10)");
	std::vector<std::string> opt_conf;
	opt_conf.push_back("field");
	opt_conf.push_back("inv-field");
	opt_conf.push_back("coeff");
	opt_conf.push_back("inv-coeff");

	bpo::variables_map vm;
	bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
	conflicting_options(vm, opt_conf);
	bpo::notify(vm);

	if (vm.count("help") || vm.size() == 0 ) {
		std::cout << all_desc << std::endl;
		return 1;
	}

	bool isField = vm.count("field") || vm.count("inv-field");
	bool isFwd = vm.count("field") || vm.count("coeff");
	bool computeInv = vm.count("compute-inverse");

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


	DisplacementFieldPointer field, field_inv;
	DisplacementFieldConstPointer input_field;

	std::vector< std::string > fnames = isFwd?fieldname:invfieldname;

	if (!isField) {
		CompositeTransformPointer tf_from_coeff = CompositeTransform::New();
		tf_from_coeff->SetPhysicalDomainInformation(ref);

		for (size_t i = 0; i < fnames.size(); i++) {
			DisplacementFieldReaderPointer fread = DisplacementFieldReaderType::New();
			fread->SetFileName( fnames[i] );
			fread->Update();
			tf_from_coeff->PushBackCoefficients(fread->GetOutput());
		}
		tf_from_coeff->Update();
		input_field = tf_from_coeff->GetDisplacementField();
	} else {
		DisplacementFieldReaderPointer fread = DisplacementFieldReaderType::New();
		fread->SetFileName( fnames[0] );
		fread->Update();
		input_field = fread->GetOutput();
	}

	const VectorType* ofb;
	VectorType* ifb;

	field = FieldType::New();
	field->SetRegions( input_field->GetLargestPossibleRegion());
	field->SetOrigin( input_field->GetOrigin());
	field->SetDirection( input_field->GetDirection() );
	field->Allocate();

	field_inv = FieldType::New();
	field_inv->SetRegions( input_field->GetLargestPossibleRegion());
	field_inv->SetOrigin( input_field->GetOrigin());
	field_inv->SetDirection( input_field->GetDirection() );
	field_inv->Allocate();

	if( isFwd ) {
		itk::ImageAlgorithm::Copy<FieldType, FieldType>(input_field, field,
				input_field->GetLargestPossibleRegion(),
				field->GetLargestPossibleRegion());
		ofb = field->GetBufferPointer();
		ifb = field_inv->GetBufferPointer();
	} else {
		itk::ImageAlgorithm::Copy<FieldType, FieldType>(input_field, field_inv,
				input_field->GetLargestPossibleRegion(),
				field_inv->GetLargestPossibleRegion());
		ofb = field_inv->GetBufferPointer();
		ifb = field->GetBufferPointer();
	}

	size_t nPix = field->GetLargestPossibleRegion().GetNumberOfPixels();
	VectorType v;
	for( size_t i = 0; i<nPix; i++ ) {
		v = *( ofb+i );
		*( ifb + i ) = -v;
	}

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
	BSplineTransformPointer tf_mesh;

	if( movingSurfaceNames.size() > 0 ){
			if(isField) {
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
			} else {
				// implement me!
			}


			tf_mesh = BSplineTransform::New();
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
}

