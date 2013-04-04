// --------------------------------------------------------------------------------------
// File:          FreeSurferTest.cxx
// Date:          Mar 22, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
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


/*
 * Test 1 - Read vtk mesh, save it as FreeSurfer binary mesh and open it with tksurfer
 * Test 2 - Read previous FreeSurfer mesh, launch registration and save result as FreeSurfer binary mesh.
 */
#ifndef DATA_DIR
#define DATA_DIR "../Data/FreeSurferTest/"
#endif

#include <itkQuadEdgeMesh.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include "MahalanobisLevelSets.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Vector<float, 2u>                               VectorPixelType;
	typedef itk::Image<VectorPixelType, 3u>                      ImageType;
	typedef MahalanobisLevelSets<ImageType>                      LevelSetsType;
	typedef LevelSetsType::ContourDeformationType                ContourDeformationType;
	typedef ContourDeformationType::Pointer                      ContourDisplacementFieldPointer;
	typedef itk::VTKPolyDataReader< ContourDeformationType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDeformationType >     WriterType;

	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( DATA_DIR ) + "new_priors_wm.vtk" );
	polyDataReader->Update();
	ContourDisplacementFieldPointer initialContour = polyDataReader->GetOutput();

	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( initialContour );
	polyDataWriter->SetFileName( std::string( DATA_DIR ) + "new_priors_wm.fsb" );

	polyDataWriter->Update();
}
