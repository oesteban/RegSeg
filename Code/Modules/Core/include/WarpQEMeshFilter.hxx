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

#ifndef WARPQEMESHFILTER_HXX_
#define WARPQEMESHFILTER_HXX_

#include "WarpQEMeshFilter.h"
#include <itkMacro.h>

namespace rstk {

template< class TInputMesh, class TOutputMesh, class TVectorField >
WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >
::WarpQEMeshFilter() {
	// Setup the number of required inputs.
	// This filter requires one mesh and one vector image as inputs.
	this->SetNumberOfRequiredInputs( 2 );
}


template< class TInputMesh, class TOutputMesh, class TVectorField >
void
WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >
::SetField(const FieldType *field ) {
	// const cast is needed because the pipeline is not const-correct.
	FieldType *input = const_cast< FieldType* >( field );
	this->itk::ProcessObject::SetNthInput( 1, input );
}


template< class TInputMesh, class TOutputMesh, class TVectorField >
const typename WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >::FieldType *
WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >
::GetField() const {
	return itkDynamicCastInDebugMode< const FieldType* >( this->itk::ProcessObject::GetInput(1) );
}

template< class TInputMesh, class TOutputMesh, class TVectorField >
void
WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >
::PrintSelf( std::ostream &os, itk::Indent indent ) const {
	Superclass::PrintSelf( os, indent );
}

/**
 * This method generates the output.
 */

template< class TInputMesh, class TOutputMesh, class TVectorField >
void
WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >
::GenerateData() {
	typedef typename TInputMesh::PointsContainer InputPointsContainer;
	typedef typename TOutputMesh::PointsContainer OutputPointsContainer;

	typedef typename TInputMesh::PointsContainerPointer InputPointsContainerPointer;
	typedef typename TOutputMesh::PointsContainerPointer OutputPointsContainerPointer;

	const InputMeshType *inputMesh = this->GetInput();
	OutputMeshPointer   outputMesh = this->GetOutput();
	FieldConstPointer        field = this->GetField();

	if ( !inputMesh ) {
		itkExceptionMacro( << "Missing Input Mesh" );
	}
	if ( !outputMesh ) {
		itkExceptionMacro( << "Missing Output Mesh" );
	}

	outputMesh->SetBufferedRegion( outputMesh->GetRequestedRegion() );

	const InputPointsContainer   *inPoints = inputMesh->GetPoints();
	OutputPointsContainerPointer outPoints = outputMesh->GetPoints();

	outPoints->Reserve( inputMesh->GetNumberOfPoints() );
	outPoints->Squeeze(); // just in case the previous mesh had allocated larger memory

	typedef typename InputMeshType::PointType InputPointType;
	typedef typename OutputMeshType::PointType OutputPointType;

	typename InputPointsContainer::ConstIterator inputPoint = inPoints->Begin();
	typename OutputPointsContainer::Iterator outputPoint = outPoints->Begin();

	size_t id = 0;
	FieldVectorType disp;
	disp.Fill(0.0);

	while (inputPoint != inPoints->End()) {
		const InputPointType p = inputPoint.Value();
		id = inputPoint.Index();
		disp = field->GetOffGridValue( id );
		outputPoint.Value() = p + disp;

		++inputPoint;
		++outputPoint;
	}

	// Create duplicate references to the rest of data on the mesh
	this->CopyInputMeshToOutputMeshPointData();
	this->CopyInputMeshToOutputMeshCells();
	this->CopyInputMeshToOutputMeshCellLinks();
	this->CopyInputMeshToOutputMeshCellData();

	// Copy Edge Cells
	typedef typename TInputMesh::CellsContainer InputCellsContainer;
	typedef typename InputCellsContainer::ConstPointer InputCellsContainerConstPointer;
	typedef typename InputCellsContainer::ConstIterator InputCellsContainerConstIterator;
	typedef typename TInputMesh::EdgeCellType InputEdgeCellType;

	InputCellsContainerConstPointer inEdgeCells = inputMesh->GetEdgeCells();

	if (inEdgeCells) {
		InputCellsContainerConstIterator ecIt = inEdgeCells->Begin();
		InputCellsContainerConstIterator ecEnd = inEdgeCells->End();

		while (ecIt != ecEnd) {
			InputEdgeCellType *pe =
					dynamic_cast<InputEdgeCellType *>(ecIt.Value());
			if (pe) {
				outputMesh->AddEdgeWithSecurePointList(pe->GetQEGeom()->GetOrigin(),
						pe->GetQEGeom()->GetDestination());
			}
			++ecIt;
		}
	}
}


}

#endif /* WARPQEMESHFILTER_HXX_ */
