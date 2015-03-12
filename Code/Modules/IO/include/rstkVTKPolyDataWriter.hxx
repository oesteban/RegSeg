// --------------------------------------------------------------------------------------
// File:          rstkVTKPolyDataWriter.hxx
// Date:          Mar 17, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
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
