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
