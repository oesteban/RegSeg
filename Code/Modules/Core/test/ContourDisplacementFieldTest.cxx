// --------------------------------------------------------------------------
// File:             ContourDisplacementFieldTest.cxx
// Date:             29/10/2012
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

#ifndef DATA_DIR
#define DATA_DIR "./"
// data source http://code.google.com/p/v3dsolver/source/browse/data/
#endif


#include "ContourDisplacementField.h"
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include <itkNormalQuadEdgeMeshFilter.h>

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef ContourDisplacementField<>             ContourDF;
	typedef typename ContourDF::Pointer            ContourDFPointer;
	typedef typename ContourDF::InternalMeshType   MeshType;
	typedef itk::VTKPolyDataReader< ContourDF >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDF >     WriterType;
	typedef itk::NormalQuadEdgeMeshFilter< ContourDF, ContourDF > NormalsFilter;

	ReaderType::Pointer polyDataReader = ReaderType::New();
	polyDataReader->SetFileName( std::string( DATA_DIR ) + "sphere.vtk" );
	polyDataReader->Update();
	ContourDFPointer mesh = polyDataReader->GetOutput();
	mesh->Initialize();

	typename ContourDF::PointsContainerPointer    points2 = mesh->GetPoints();
	typename ContourDF::PointsContainerIterator   p_it2 = points2->Begin();
	typename ContourDF::PointDataContainerPointer container2 = mesh->GetPointData();
	typename ContourDF::PointDataContainerIterator d_it2 = container2->Begin();

	  std::cout << "Index2 * Point2 * Data2" << std::endl;

	  while( p_it2 != points2->End() )
	    {
	    std::cout << p_it2.Index() << " * ";
	    std::cout << p_it2.Value() << " * ";
	    std::cout << d_it2.Value() << std::endl;

	    ++p_it2;
	    ++d_it2;
	    }

		WriterType::Pointer polyDataWriter = WriterType::New();
		polyDataWriter->SetInput( mesh );
		polyDataWriter->SetFileName( std::string( DATA_DIR ) + "sphere-field.vtp" );
		polyDataWriter->Update();

	/*
	typename NormalsFilter::Pointer norm = NormalsFilter::New();
	norm->SetInput( mesh );
	norm->Update();

	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( norm->GetOutput() );
	polyDataWriter->SetFileName( std::string( DATA_DIR ) + "sphere-normals.vtk" );
	polyDataWriter->Update();

	typename NormalsFilter::OutputMeshType::Pointer output = norm->GetOutput( );

	typename NormalsFilter::OutputMeshType::PointsContainerPointer    points = output->GetPoints();
	typename NormalsFilter::OutputMeshType::PointsContainerIterator   p_it = points->Begin();

	typename NormalsFilter::OutputMeshType::PointDataContainerPointer container = output->GetPointData();
	typename NormalsFilter::OutputMeshType::PointDataContainerIterator d_it = container->Begin();

	  std::cout << "Index * Point * Normal" << std::endl;

	  while( p_it != points->End() )
	    {
	    std::cout << p_it.Index() << " * ";
	    std::cout << p_it.Value() << " * ";
	    std::cout << d_it.Value() << std::endl;

	    ++p_it;
	    ++d_it;
	    }
	    */

/*
  	typedef typename ContourDF::PointDataContainer::Iterator PointsIterator;
	PointsIterator  pointIterator = mesh->GetPoints()->Begin();
	PointsIterator end = mesh->GetPoints()->End();
	  int nPoints = 0;
	  while( pointIterator != end ) {
		  p = pointIterator.;

	    if( p != pts[ nPoints ] )
	      {
	      std::cout << "Point N. " << nPoints << " differs." << std::endl;
	      return EXIT_FAILURE;
	      }

	    pointIterator++;
	    nPoints++;
	    }
*/
	//ContourDFPointer p = ContourDF::New();
	//p->SetReferenceMesh( polyDataReader->GetOutput() );
}
