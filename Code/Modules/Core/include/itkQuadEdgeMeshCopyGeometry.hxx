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
