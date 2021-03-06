// This file is part of RegSeg
//
// Copyright 2014-2017, Oscar Esteban <code@oscaresteban.es>
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

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
			("write-field,w", bpo::bool_switch(), "write output field")
			("mask,m", bpo::value< std::string >(&maskfile), "mask file" )
			("mask-inputs", bpo::bool_switch(), "use deformed mask to filter input files")
			("field,F", bpo::value < std::vector< std::string > >(&fieldname), "forward displacement field" )
			("inv-field,R", bpo::value < std::vector< std::string > >(&invfieldname), "backward displacement field" )
			("coeff,C", bpo::value < std::vector< std::string > >(&fieldname)->multitoken(), "forward displacement field" )
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
	DisplacementFieldPointer input_field;

	std::vector< std::string > fnames = isFwd?fieldname:invfieldname;
	std::vector< CoefficientsImageArray > allcoeff;

	if (!isField) { // Read coefficients
		for (size_t i = 0; i < fnames.size(); i++) {
			CoefficientsImageArray coeffarr;
			VectorFieldReaderPointer fread = VectorFieldReaderType::New();
			fread->SetFileName( fnames[i] );
			fread->Update();
			FakeFieldType::Pointer f = fread->GetOutput();
			size_t tnp = f->GetLargestPossibleRegion().GetNumberOfPixels();

			CoefficientsImageType::SizeType size;
			CoefficientsImageType::SpacingType sp;
			CoefficientsImageType::PointType orig;

			for (size_t dim = 0; dim < DIMENSION; dim++) {
				size[dim] = f->GetLargestPossibleRegion().GetSize()[dim];
				sp[dim] = f->GetSpacing()[dim];
				orig[dim] = - 0.5 * size[dim] * sp[dim];
			}

			float* dbuff[DIMENSION];
			for (size_t dim = 0; dim < DIMENSION; dim++) {
				coeffarr[dim] = CoefficientsImageType::New();
				coeffarr[dim]->SetRegions(size);
				coeffarr[dim]->SetSpacing(sp);
				coeffarr[dim]->SetOrigin(orig);
				coeffarr[dim]->Allocate();
				dbuff[dim] = coeffarr[dim]->GetBufferPointer();
			}
			size_t np = coeffarr[0]->GetLargestPossibleRegion().GetNumberOfPixels();

			float* sbuff = f->GetBufferPointer();
			float svect;
			size_t pixid, dim;
			for( size_t pix = 0; pix < tnp; pix++) {
				svect = *(sbuff + pix);
				pixid = pix % np;
				dim = pix / np;
				*( dbuff[dim] + pixid ) = svect;
			}
			allcoeff.push_back(coeffarr);
		}
	}

	if(!isField) {
		CompositeTransformPointer tf_from_coeff = CompositeTransform::New();
		tf_from_coeff->SetOutputReference(ref);
		for (size_t i = 0; i < allcoeff.size(); i++) {
			tf_from_coeff->PushBackCoefficients(allcoeff[i]);
		}
		tf_from_coeff->Interpolate();
		input_field = tf_from_coeff->GetDisplacementField();
	} else {
		VectorFieldReaderPointer fread = VectorFieldReaderType::New();
		fread->SetFileName( fnames[0] );
		fread->Update();
		FakeFieldType::Pointer f = fread->GetOutput();
		size_t tnp = f->GetLargestPossibleRegion().GetNumberOfPixels();

		FieldType::SizeType size;
		FieldType::SpacingType sp;
		FieldType::PointType orig;
		FieldType::DirectionType dir;
		dir.Fill(0.0);

		for (size_t dim = 0; dim < DIMENSION; dim++) {
			size[dim] = f->GetLargestPossibleRegion().GetSize()[dim];
			sp[dim] = f->GetSpacing()[dim];
			orig[dim] = f->GetOrigin()[dim];
			for (size_t j = 0; j < DIMENSION; j++) {
				dir[dim][j] = f->GetDirection()[dim][j];
			}
		}


		VectorType v;
		v.Fill(0.0);
		input_field = FieldType::New();
		input_field->SetRegions(size);
		input_field->SetSpacing(sp);
		input_field->SetOrigin(orig);
		input_field->SetDirection(dir);
		input_field->Allocate();
		input_field->FillBuffer(v);

		size_t np = input_field->GetLargestPossibleRegion().GetNumberOfPixels();

		float* sbuff = f->GetBufferPointer();
		float svect;
		VectorType* dbuff = input_field->GetBufferPointer();
		size_t pixid, dim;
		for( size_t pix = 0; pix < tnp; pix++) {
			svect = *(sbuff + pix);
			pixid = pix % np;
			dim = pix / np;
			v = *(dbuff + pixid);
			v[dim] = svect;
			*( dbuff + pixid ) = v;
		}
	}

	const VectorType* ofb;
	VectorType* ifb;

	field = FieldType::New();
	field->SetRegions( input_field->GetLargestPossibleRegion());
	field->SetOrigin( input_field->GetOrigin());
	field->SetSpacing( input_field->GetSpacing() );
	field->SetDirection( input_field->GetDirection() );
	field->Allocate();

	field_inv = FieldType::New();
	field_inv->SetRegions( input_field->GetLargestPossibleRegion());
	field_inv->SetSpacing( input_field->GetSpacing() );
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

	if (vm.count("write-field") && vm["write-field"].as<bool>() ) {
		FieldWriterPointer w = FieldWriter::New();
		w->SetInput(field);
		w->SetFileName((outPrefix + "_field.nii.gz").c_str());
		w->Update();
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
	BaseTransformPointer tf_mesh;

	if( movingSurfaceNames.size() > 0 ){
			if(isField) {
				if ( vm.count("grid-size") ){
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
				BSplineTransformPointer tfm = BSplineTransform::New();
				tfm->SetControlGridSize(size);
				tfm->SetDomainExtent( field );
				tfm->SetDisplacementField( field );
				tfm->ComputeCoefficients();
				tf_mesh = itkDynamicCastInDebugMode< BaseTransform* >(tfm.GetPointer());
			} else {
				CompositeTransformPointer tfm = CompositeTransform::New();

				for (size_t i = 0; i < allcoeff.size(); i++) {
					tfm->PushBackCoefficients(allcoeff[i]);
				}
				tf_mesh = itkDynamicCastInDebugMode< BaseTransform* >(tfm.GetPointer());
			}

			PointsList points;
			for( size_t i = 0; i<movingSurfaceNames.size(); i++){
				MeshReaderPointer r = MeshReaderType::New();
				r->SetFileName( movingSurfaceNames[i] );
				r->Update();

				MeshPointer cur_mesh = r->GetOutput();
				PointsIterator p_it = cur_mesh->GetPoints()->Begin();
				PointsIterator p_end = cur_mesh->GetPoints()->End();

				MeshPointType p;
				while ( p_it!=p_end ) {
					points.push_back(p_it.Value());
					++p_it;
				}
			}
			tf_mesh->SetOutputPoints(points);
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
					p_it.Value() += tf_mesh->GetPointValue(pointId);

					++p_it;
					pointId++;
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

