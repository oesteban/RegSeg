// --------------------------------------------------------------------------
// File:             MahalanobisLevelSetsTest.cxx
// Date:             30/10/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWEReg
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "./"
// data source http://code.google.com/p/v3dsolver/source/browse/data/
#endif

#include <itkVector.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include "MahalanobisLevelSets.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::QuadEdgeMesh< itk::Vector< double, 3u >, 3u > ContourDisplacementFieldType;
	typedef typename ContourDisplacementFieldType::Pointer     ContourDisplacementFieldPointer;

	typedef itk::VTKPolyDataReader< ContourDisplacementFieldType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDisplacementFieldType >     WriterType;


	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( TEST_DATA_DIR ) + "sphere.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer mesh = polyDataReader->GetOutput();
	mesh->Initialize();



		WriterType::Pointer polyDataWriter = WriterType::New();
		polyDataWriter->SetInput( mesh );
		polyDataWriter->SetFileName( std::string( TEST_DATA_DIR ) + "sphere-field.vtp" );
		polyDataWriter->Update();

}
