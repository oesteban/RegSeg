/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadEdgeMeshVectorDataVTKPolyDataWriter.txx,v $
  Language:  C++
  Date:      $Date: 2009-01-02 18:43:05 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkQuadEdgeMeshVectorDataVTKPolyDataWriter_txx
#define __itkQuadEdgeMeshVectorDataVTKPolyDataWriter_txx

#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include <itkMeasurementVectorTraits.h>


namespace itk
{

// 
// Constructor
// 
template< class TMesh >
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh> 
::QuadEdgeMeshVectorDataVTKPolyDataWriter() 
{
  m_CellDataName = "celldatascalar";
  m_PointDataName = "pointdatavecs";
}

// 
// Destructor
// 
template< class TMesh >
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh> 
::~QuadEdgeMeshVectorDataVTKPolyDataWriter()
{
}

template< class TMesh >
void
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh> 
::GenerateData()
{
  this->Superclass::GenerateData();
  this->WriteCellData();
  this->WritePointData();
}

template< class TMesh >
void
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh> 
::WriteCellData()
{
  CellDataContainerConstPointer celldata = this->m_Input->GetCellData();

  if( celldata )
    {
    if( celldata->Size() != 0 )
      {
      std::ofstream outputFile( this->m_FileName.c_str(), std::ios_base::app );

      outputFile <<"CELL_DATA " <<this->m_Input->GetNumberOfFaces() <<std::endl;
      outputFile <<"SCALARS ";
      outputFile <<m_CellDataName <<" double 1" <<std::endl;
      outputFile <<"LOOKUP_TABLE default" <<std::endl;

      CellsContainerConstPointer cells = this->m_Input->GetCells();
      CellsContainerConstIterator it = cells->Begin();

      CellDataContainerConstIterator c_it = celldata->Begin();

      while( c_it != celldata->End() )
        {
        CellType* cellPointer = it.Value();
        outputFile <<c_it.Value();
        ++c_it;
        ++it;
        }
      outputFile <<std::endl;
      outputFile.close();
      }
    }
}

template< class TMesh >
void
QuadEdgeMeshVectorDataVTKPolyDataWriter<TMesh> 
::WritePointData()
{
  PointDataContainerConstPointer pointdata = this->m_Input->GetPointData();

  if( pointdata )
    {
    std::ofstream outputFile( this->m_FileName.c_str(), std::ios_base::app );

    outputFile <<"POINT_DATA " <<this->m_Input->GetNumberOfPoints() <<std::endl;
    outputFile <<"vectors ";

    PointDataContainerIterator c_it = pointdata->Begin();
    typedef typename MeshType::PointDataContainer::Element VectorType;
    typedef NumericTraits<VectorType> ElementTraits;
    const unsigned int vectorSize =  ElementTraits::GetLength( c_it.Value() );
    outputFile << m_PointDataName << " double" << std::endl;

    VectorType v;
    while(  c_it != pointdata->End() )
      {
      for( unsigned int j = 0; j < vectorSize; j++ )
        {
    	  v = c_it.Value();
    	  outputFile << static_cast< double >( v[j] );

    	  if( (j + 1) % vectorSize == 0 ) {
    		  outputFile << std::endl;
    	  } else {
    		  outputFile << " ";
    	  }
        }

      ++c_it;
      }

    outputFile << std::endl;
    outputFile.close();
    }
}

}

#endif
