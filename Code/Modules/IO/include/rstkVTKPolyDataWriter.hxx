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

#ifndef RSTKVTKPOLYDATAWRITER_HXX_
#define RSTKVTKPOLYDATAWRITER_HXX_

#include "rstkVTKPolyDataWriter.h"
#include "itkCellInterface.h"
#include <fstream>

namespace rstk {

template< class TInputMesh >
void
VTKPolyDataWriter< TInputMesh >
::GenerateData()
{
  if ( this->m_FileName == "" )
    {
    itkExceptionMacro("No FileName");
    return;
    }

  //
  // Write to output file
  //
  std::ofstream outputFile( this->m_FileName.c_str() );

  if ( !outputFile.is_open() )
    {
    itkExceptionMacro("Unable to open file\n"
                      "outputFilename= " << this->m_FileName);
    return;
    }

  outputFile.imbue( std::locale::classic() );
  outputFile << "# vtk DataFile Version 1.0" << std::endl;
  outputFile << "vtk output" << std::endl;
  outputFile << "ASCII" << std::endl;
  outputFile << "DATASET POLYDATA" << std::endl;

  // POINTS go first

  unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();
  outputFile << "POINTS " << numberOfPoints << " float" << std::endl;

  const PointsContainer *points = this->m_Input->GetPoints();

  std::map< PointIdentifier, PointIdentifier > IdMap;
  PointIdentifier                              k = 0;

  if ( points )
    {
    PointIterator pointIterator = points->Begin();
    PointIterator pointEnd = points->End();

    while ( pointIterator != pointEnd )
      {
      PointType point = pointIterator.Value();

      outputFile << point[0] << " " << point[1];

      if ( TInputMesh::PointDimension > 2 )
        {
        outputFile << " " << point[2];
        }
      else
        {
        outputFile << " " << "0.0";
        }

      outputFile << std::endl;

      IdMap[pointIterator.Index()] = k++;
      pointIterator++;
      }
    }

  unsigned int numberOfVertices = 0;
  unsigned int numberOfEdges = 0;
  unsigned int numberOfPolygons = 0;

  const CellsContainer *cells = this->m_Input->GetCells();

  if ( cells )
    {
    CellIterator cellIterator = cells->Begin();
    CellIterator cellEnd = cells->End();

    while ( cellIterator != cellEnd )
      {
      switch ( cellIterator.Value()->GetType() )
        {
        case 0: //VERTEX_CELL:
          numberOfVertices++;
          break;
        case 1: //LINE_CELL:
        case 7: //QUADRATIC_EDGE_CELL:
          numberOfEdges++;
          break;
        case 2: //TRIANGLE_CELL:
        case 3: //QUADRILATERAL_CELL:
        case 4: //POLYGON_CELL:
        case 8: //QUADRATIC_TRIANGLE_CELL:
          numberOfPolygons++;
          break;
        default:
          std::cerr << "Unhandled cell (volumic?)." << std::endl;
          break;
        }
      cellIterator++;
      }

    // VERTICES should go here
    if ( numberOfVertices )
        {}

    // LINES
    if ( numberOfEdges )
      {
      outputFile << "LINES " << numberOfEdges << " " << 3 * numberOfEdges;
      outputFile << std::endl;

      cellIterator = cells->Begin();
      while ( cellIterator != cellEnd )
        {
        CellType *cellPointer = cellIterator.Value();
        switch ( cellIterator.Value()->GetType() )
          {
          case 1: //LINE_CELL:
          case 7: //QUADRATIC_EDGE_CELL:
            {
            PointIdIterator pointIdIterator = cellPointer->PointIdsBegin();
            PointIdIterator pointIdEnd = cellPointer->PointIdsEnd();
            outputFile << cellPointer->GetNumberOfPoints();
            while ( pointIdIterator != pointIdEnd )
              {
              outputFile << " " << IdMap[*pointIdIterator];
              pointIdIterator++;
              }
            outputFile << std::endl;
            break;
            }
          default:
            break;
          }
        cellIterator++;
        }
      }

    // POLYGONS
    if ( numberOfPolygons )
      {
      // This could be optimized but at least now any polygonal
      // mesh can be saved.
      cellIterator = cells->Begin();

      PointIdentifier totalNumberOfPointsInPolygons = itk::NumericTraits<PointIdentifier>::Zero;
      while ( cellIterator != cells->End() )
        {
        CellType *cellPointer = cellIterator.Value();
        if ( cellPointer->GetType() != CellType::VERTEX_CELL
             && cellPointer->GetType() != CellType::LINE_CELL )
          {
          totalNumberOfPointsInPolygons += cellPointer->GetNumberOfPoints();
          }
        cellIterator++;
        }
      outputFile << "POLYGONS " << numberOfPolygons << " ";
      outputFile << totalNumberOfPointsInPolygons + numberOfPolygons; // FIXME: Is this right ?
      outputFile << std::endl;

      cellIterator = cells->Begin();
      while ( cellIterator != cellEnd )
        {
        CellType *cellPointer = cellIterator.Value();
        switch ( cellIterator.Value()->GetType() )
          {
          case 2: //TRIANGLE_CELL:
          case 3: //QUADRILATERAL_CELL:
          case 4: //POLYGON_CELL:
          case 8: //QUADRATIC_TRIANGLE_CELL:
            {
            PointIdIterator pointIdIterator = cellPointer->PointIdsBegin();
            PointIdIterator pointIdEnd = cellPointer->PointIdsEnd();
            outputFile << cellPointer->GetNumberOfPoints();
            while ( pointIdIterator != pointIdEnd )
              {
              outputFile << " " << IdMap[*pointIdIterator];
              pointIdIterator++;
              }
            outputFile << std::endl;
            break;
            }
          default:
            break;
          }
        cellIterator++;
        }
      }
    }

  // TRIANGLE_STRIP should go here
  // except that ... there is no such thing in ITK ...

  outputFile.close();
}

}  // namespace rstk

#endif /* RSTKVTKPOLYDATAWRITER_HXX_ */
