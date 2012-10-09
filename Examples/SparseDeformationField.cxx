// --------------------------------------------------------------------------------------
// File:          SparseDeformationField.cxx
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

#include <itkMesh.h>
#include <itkVector.h>
#include "SparseDeformationField.h"


typedef float       PixelType;
const unsigned int Dimension = 3;
typedef itk::Vector< PixelType, 3 >                 VectorType;
typedef itk::Mesh< VectorType, Dimension >      SurfaceType;
typedef rstk::SparseDeformationField<SurfaceType, double> FieldType;

int main(int argc, char *argv[]) {
/*	unsigned int npoints = 4;
	// Create PointSet
	SurfaceType::Pointer svf = SurfaceType::New();
	svf->SetRequestedRegion( npoints );

	// Set some vectors
	SurfaceType::PointType p1;
	p1.Fill(10);
	svf->SetPoint( 0, p1 );

	SurfaceType::PointType p2;
	p2.Fill( 90.0 );
	svf->SetPoint( 1, p2 );

	SurfaceType::PointType p3;
	p3[0] = 10;
	p3[1] = 90;
	p3[2] = 90;
	svf->SetPoint( 2, p3 );

	SurfaceType::PointType p4;
	p4[0] = 90;
	p4[1] = 90;
	p4[2] = 10;
	svf->SetPoint( 3, p4 );

	// Deformation field
	FieldType::Pointer f = FieldType::New();
	f->SetReferenceSurface(svf);

	VectorType p;
	p.Fill( 1.0 );
	for ( size_t i = 0; i<npoints; i++ ) {
		f->SetK( i, p );
	}

	SurfaceType::ConstPointer tf = f->GetTransformedSurface();


	for ( size_t i = 0; i<npoints; i++ ) {
		tf->GetPoint( i, &p1 );

		std::cout << p1 << std::endl;
	}

	// Write result
	//Writer::Pointer w = Writer::New();
	//w->SetInput( res->GetOutput() );
	//w->SetFileName( "TransformedMesh.vtk");
	//w->Update();
*/
	return EXIT_SUCCESS;
}
