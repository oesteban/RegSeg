// --------------------------------------------------------------------------------------
// File:          GaussianModelGradientComputation.cxx
// Date:          Apr 5, 2012
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
//   * Neither the names of the copyright holders, nor the names of
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
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
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//


#include <itkPointSet.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkMesh.h>
#include <itkMeshFileWriter.h>
#include <itkBinaryMask3DMeshSource.h>

#include "SparseDeformationField.h"
#include "GaussianModelDataEnergy.h"

typedef float       PixelType;
const unsigned int Dimension = 3;
typedef itk::Vector< PixelType, 3 >                 VectorType;
typedef itk::Image< unsigned char, 3 >              ImageType;
typedef itk::Mesh< unsigned char, 3 >               SurfaceType;
typedef itk::MeshFileWriter< SurfaceType >          SurfaceWriterType;


typedef itk::BinaryMask3DMeshSource< ImageType, SurfaceType >   Image2MeshType;


typedef rstk::GaussianModelDataEnergy< ImageType, SurfaceType >      ModelType;
typedef rstk::SparseDeformationField< SurfaceType, double >          FieldType;


int main(int argc, char *argv[]) {
	// Generate a reference image
	ImageType::Pointer ref = ImageType::New();
	ImageType::SizeType ref_size;
	ref_size.Fill( 100 );
	ref->SetRegions( ref_size );
	ref->SetSpacing(1.0);
	ref->Allocate();
	ref->FillBuffer( 0 );
	ref->Update();

	ImageType::IndexType idx;
	for( size_t i = 24; i < 75; i++ ) {
		idx[0] = i;
		for( size_t j = 24; j < 75; j++ ) {
			idx[1] = j;
			for( size_t k = 24; k < 75; k++ ) {
				idx[2] = k;

				ref->SetPixel( idx, 1 );
			}
		}
	}

	typedef itk::ImageFileWriter< ImageType > Writer;
	Writer::Pointer w = Writer::New();
	w->SetFileName( "referenceImage.vtk" );
	w->SetInput( ref );
	w->Update();

	// Generate an initial surface
	Image2MeshType::Pointer i2m = Image2MeshType::New();
	i2m->SetInput( ref );
	i2m->SetRegionOfInterest( ref->GetLargestPossibleRegion() );
	i2m->Update();

	SurfaceWriterType::Pointer m_w = SurfaceWriterType::New();
	m_w->SetInput( i2m->GetOutput() );
	m_w->SetFileName( "movingMesh.vtk" );
	m_w->Update();

	// Generate a sparse deformation field for every point in the surface
	FieldType::Pointer field = FieldType::New();
	field->SetReferenceSurface( i2m->GetOutput() );

	VectorType displacement; displacement.Fill( 5.0 );
	displacement[2] = 3.0;

	field->SetK( 0, displacement );

	m_w->SetInput( field->GetTransformedSurface() );
	m_w->SetFileName( "movingMeshDeformed.vtk" );
	m_w->Update();


	// Set-up

	// Output gradient

	return EXIT_SUCCESS;
}
