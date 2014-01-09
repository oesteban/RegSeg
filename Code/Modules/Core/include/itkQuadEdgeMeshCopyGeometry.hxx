// --------------------------------------------------------------------------------------
// File:          itkQuadEdgeMeshCopyGeometry.hxx
// Date:          Jan 9, 2014
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

#ifndef ITKQUADEDGEMESHCOPYGEOMERTY_HXX_
#define ITKQUADEDGEMESHCOPYGEOMERTY_HXX_

#include "itkQuadEdgeMeshCopyGeometry.h"

namespace itk
{
// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::QuadEdgeMeshCopyGeometry()
{
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNumberOfRequiredOutputs(1);

  this->Superclass::SetNthOutput( 0, OutputMeshType::New() );
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyMesh()
{

	  const InputMeshType *in = this->GetInput();
	  OutputMeshType *     out = this->GetOutput();
	  this->Copy(in, out);
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyMeshGeometry()
{
  this->CopyMeshPoints();
  this->CopyMeshEdgeCells();
  this->CopyMeshCells();
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyMeshPoints()
{
  const InputMeshType *in = this->GetInput();
  OutputMeshType *     out = this->GetOutput();

  CopyPoints(in, out);
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyMeshEdgeCells()
{
  const InputMeshType *in = this->GetInput();
  OutputMeshType *     out = this->GetOutput();

  CopyEdgeCells(in, out);
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyMeshCells()
{
  const InputMeshType *in = this->GetInput();
  OutputMeshType *     out = this->GetOutput();

  CopyCells(in, out);
}

//
// Helper functions that copy selected pieces of a Mesh.
// These functions should be templated here in order to
// facilitate their reuse in multiple scenarios.
//
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::Copy(const TInputMesh *in, TOutputMesh *out)
{
  CopyPoints(in, out);
  CopyEdgeCells(in, out);
  CopyCells(in, out);

  typedef typename TInputMesh::PointDataContainer        InputPointDataContainer;
  typedef typename TOutputMesh::PointDataContainer       OutputPointDataContainer;
  typedef typename InputPointDataContainer::ConstPointer InputPointDataContainerConstPointer;
  typedef typename OutputPointDataContainer::Pointer     OutputPointDataContainerPointer;

  InputPointDataContainerConstPointer inputPointData = in->GetPointData();
  if ( inputPointData.IsNotNull() )
    {
    // There is nothing to copy
	  OutputPointDataContainerPointer outputPointData = OutputPointDataContainer::New();
	  outputPointData->Reserve( inputPointData->Size() );

	  // Copy point data
	  typedef typename InputPointDataContainer::ConstIterator InputPointDataContainerConstIterator;
	  InputPointDataContainerConstIterator inIt = inputPointData->Begin();
	  while ( inIt != inputPointData->End() )
	    {
	    outputPointData->SetElement( inIt.Index(), 0.0 );
	    inIt++;
	    }

	  out->SetPointData(outputPointData);
    }
  else {
	  size_t dataSize = out->ComputeNumberOfPoints();
	  OutputPointDataContainerPointer outputPointData = OutputPointDataContainer::New();
	  outputPointData->Reserve( dataSize );

	  for( size_t i = 0; i<dataSize; i++ ){
		  outputPointData->SetElement( i , 0.0 );
	  }
	  out->SetPointData( outputPointData );

  }

}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyCells(const TInputMesh *in, TOutputMesh *out)
{
  // Copy cells
  typedef typename TInputMesh::CellsContainer         InputCellsContainer;
  typedef typename InputCellsContainer::ConstPointer  InputCellsContainerConstPointer;
  typedef typename InputCellsContainer::ConstIterator InputCellsContainerConstIterator;
  typedef typename TInputMesh::PolygonCellType        InputPolygonCellType;
  typedef typename TInputMesh::PointIdList            InputPointIdList;
  typedef typename TInputMesh::CellTraits             InputCellTraits;
  typedef typename InputCellTraits::PointIdInternalIterator
  InputPointsIdInternalIterator;

  out->SetCellsAllocationMethod(TOutputMesh::CellsAllocatedDynamicallyCellByCell);

  InputCellsContainerConstPointer inCells = in->GetCells();

  if ( inCells )
    {
    InputCellsContainerConstIterator cIt = inCells->Begin();
    InputCellsContainerConstIterator cEnd = inCells->End();
    while ( cIt != cEnd )
      {
      InputPolygonCellType *pe = dynamic_cast< InputPolygonCellType * >( cIt.Value() );
      if ( pe )
        {
        InputPointIdList              points;
        InputPointsIdInternalIterator pIt   = pe->InternalPointIdsBegin();
        InputPointsIdInternalIterator pEnd  = pe->InternalPointIdsEnd();

        while ( pIt != pEnd )
          {
          points.push_back( ( *pIt ) );
          ++pIt;
          }
        out->AddFaceWithSecurePointList(points, false);
        }
      ++cIt;
      }
    }
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyEdgeCells(const TInputMesh *in, TOutputMesh *out)
{
  // Copy Edge Cells
  typedef typename TInputMesh::CellsContainer         InputCellsContainer;
  typedef typename InputCellsContainer::ConstPointer  InputCellsContainerConstPointer;
  typedef typename InputCellsContainer::ConstIterator InputCellsContainerConstIterator;
  typedef typename TInputMesh::EdgeCellType           InputEdgeCellType;

  InputCellsContainerConstPointer inEdgeCells = in->GetEdgeCells();

  if ( inEdgeCells )
    {
    InputCellsContainerConstIterator ecIt   = inEdgeCells->Begin();
    InputCellsContainerConstIterator ecEnd  = inEdgeCells->End();

    while ( ecIt != ecEnd )
      {
      InputEdgeCellType *pe = dynamic_cast< InputEdgeCellType * >( ecIt.Value() );
      if ( pe )
        {
        out->AddEdgeWithSecurePointList( pe->GetQEGeom()->GetOrigin(),
                                         pe->GetQEGeom()->GetDestination() );
        }
      ++ecIt;
      }
    }
}

// ---------------------------------------------------------------------
template< class TInputMesh, class TOutputMesh >
void
QuadEdgeMeshCopyGeometry< TInputMesh, TOutputMesh >
::CopyPoints(const TInputMesh *in, TOutputMesh *out)
{
  // Copy points
  typedef typename TInputMesh::PointsContainerConstPointer   InputPointsContainerConstPointer;
  typedef typename TInputMesh::PointsContainerConstIterator  InputPointsContainerConstIterator;

  typedef typename TOutputMesh::PointsContainer             OutputPointsContainer;
  typedef typename TOutputMesh::PointsContainerPointer      OutputPointsContainerPointer;
  typedef typename TOutputMesh::PointType                   OutputPointType;

  InputPointsContainerConstPointer inPoints = in->GetPoints();

  if ( inPoints )
    {
    InputPointsContainerConstIterator inIt  = inPoints->Begin();
    InputPointsContainerConstIterator inEnd = inPoints->End();

    OutputPointsContainerPointer      oPoints = out->GetPoints();
    if( oPoints.IsNull() )
      {
      oPoints = OutputPointsContainer::New();
      out->SetPoints( oPoints );
      }
    OutputPointType                   pOut;

    while ( inIt != inEnd )
      {
      pOut.CastFrom( inIt.Value() );
      oPoints->InsertElement(inIt.Index(), pOut);
      ++inIt;
      }
    }
}

} // end namespace itk

#endif /* ITKQUADEDGEMESHCOPYGEOMERTY_HXX_ */
