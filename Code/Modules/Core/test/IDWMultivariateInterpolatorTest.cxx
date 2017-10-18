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

#include "SparseToDenseFieldResampleFilter.h"
#include "IDWMultivariateInterpolator.h"


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMesh.h>
#include <itkMeshFileReader.h>
#include "DisplacementFieldFileWriter.h"


namespace bpo = boost::program_options;

typedef float       PixelType;
const unsigned int Dimension = 3;
typedef itk::Vector< PixelType, 3 >                 PointType;
typedef itk::Image< PointType, Dimension >         DenseVectorFieldType;
typedef itk::PointSet< PointType, Dimension >      SparseVectorFieldType;
typedef rstk::SparseToDenseFieldResampleFilter<SparseVectorFieldType, DenseVectorFieldType>  ResamplerType;
//typedef itk::ImageFileWriter< DenseVectorFieldType > Writer;
typedef rstk::DisplacementFieldFileWriter<DenseVectorFieldType> Writer;

int main(int argc, char *argv[]) {
	// Create PointSet
	SparseVectorFieldType::Pointer svf = SparseVectorFieldType::New();
	svf->SetRequestedRegion( 3 );


	// Set some vectors
	PointType p1;
	PointType v1;
	v1.Fill( 5.0 );
	p1[0] = 10;
	p1[1] = 10;
	p1[2] = 10;
	svf->SetPoint( 0, p1 );
	svf->SetPointData( 0, v1 );

	PointType p2;
	p2[0] = 90;
	p2[1] = 90;
	p2[2] = 90;
	PointType v2;
	v2[0] = -5.0;
	v2[1] = -15.0;
	v2[2] = 10.0;

	svf->SetPoint( 1, p2 );
	svf->SetPointData( 1, v2 );

	PointType p3;
	p3[0] = 10;
	p3[1] = 90;
	p3[2] = 70;
	PointType v3;
	v3[0] = 2.0;
	v3[1] = -1.0;
	v3[2] = -5.0;

	svf->SetPoint( 2, p3 );
	svf->SetPointData( 2, v3 );

	// Resample
	ResamplerType::OutputImageSizeType size;
	size.Fill( 50 );
	ResamplerType::Pointer res = ResamplerType::New();
	res->SetInput( svf );
	res->AddControlPoints( svf );
	res->SetFieldSize( size );
	res->SetFieldSpacing( 2.0 );
	res->Update();
	DenseVectorFieldType::Pointer df = res->GetOutput();

	// Write result
	Writer::Pointer p = Writer::New();
	p->SetFileName( "IDWMultivariateInterpolatorTest.nii.gz" );
	p->SetInput( df );
	p->Update();

	return EXIT_SUCCESS;
}
