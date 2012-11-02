// --------------------------------------------------------------------------------------
// File:          IDWMultivariateInterpolatorTest.cxx
// Date:          Mar 29, 2012
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       0.1
// License:       BSD
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Biomedical Image Technology, UPM (BIT-UPM)
// and Signal Processing Lab 5, EPFL (LTS5-EPFL)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//   * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names
// of its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//


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
#include <itkImageToVectorImageFilter.h>
#include "DisplacementFieldFileWriter.h"


namespace bpo = boost::program_options;

typedef float       PixelType;
const unsigned int Dimension = 3;
typedef itk::Vector< PixelType, 3 >                 VectorType;
typedef itk::Image< VectorType, Dimension >         DenseVectorFieldType;
typedef itk::PointSet< VectorType, Dimension >      SparseVectorFieldType;
typedef rstk::SparseToDenseFieldResampleFilter<SparseVectorFieldType, DenseVectorFieldType>  ResamplerType;
//typedef itk::ImageFileWriter< DenseVectorFieldType > Writer;
typedef rstk::DisplacementFieldFileWriter<DenseVectorFieldType> Writer;

int main(int argc, char *argv[]) {
	// Create PointSet
	SparseVectorFieldType::Pointer svf = SparseVectorFieldType::New();
	svf->SetRequestedRegion( 3 );


	// Set some vectors
	SparseVectorFieldType::PointType p1;
	VectorType v1;
	v1.Fill( 5.0 );
	p1[0] = 10;
	p1[1] = 10;
	p1[2] = 10;
	svf->SetPoint( 0, p1 );
	svf->SetPointData( 0, v1 );

	SparseVectorFieldType::PointType p2;
	p2[0] = 90;
	p2[1] = 90;
	p2[2] = 90;
	VectorType v2;
	v2[0] = -5.0;
	v2[1] = -15.0;
	v2[2] = 10.0;

	svf->SetPoint( 1, p2 );
	svf->SetPointData( 1, v2 );

	SparseVectorFieldType::PointType p3;
	p3[0] = 10;
	p3[1] = 90;
	p3[2] = 70;
	VectorType v3;
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
	res->SetOutputSize( size );
	res->SetOutputSpacing( 2.0 );
	res->Update();
	DenseVectorFieldType::Pointer df = res->GetOutput();

	// Write result
	Writer::Pointer p = Writer::New();
	p->SetFileName( "IDWMultivariateInterpolatorTest.nii.gz" );
	p->SetInput( df );
	p->Update();

	return EXIT_SUCCESS;
}
