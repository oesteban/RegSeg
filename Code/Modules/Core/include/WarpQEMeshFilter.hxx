// --------------------------------------------------------------------------------------
// File:          WarpQEMeshFilter.hxx
// Date:          Jan 14, 2014
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
	this->ProcessObject::SetNthInput( 1, input );
}


template< class TInputMesh, class TOutputMesh, class TVectorField >
const typename WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >::FieldType *
WarpQEMeshFilter< TInputMesh, TOutputMesh, TVectorField >
::GetField() {
	return itkDynamicCastInDebugMode< const FieldType* >( this->ProcessObject::GetInput(1) );
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
	typedef typename FieldType::IndexType IndexType;
	IndexType index;

	const unsigned int Dimension = field->GetImageDimension();

	const InputPointType p;
	VectorType disp;
	while (inputPoint != inPoints->End()) {
		p = inputPoint.Value();
		disp = field->GetOffGridValue( inputPoint.Index() );
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

	InputCellsContainerConstPointer inEdgeCells = in->GetEdgeCells();

	if (inEdgeCells) {
		InputCellsContainerConstIterator ecIt = inEdgeCells->Begin();
		InputCellsContainerConstIterator ecEnd = inEdgeCells->End();

		while (ecIt != ecEnd) {
			InputEdgeCellType *pe =
					dynamic_cast<InputEdgeCellType *>(ecIt.Value());
			if (pe) {
				out->AddEdgeWithSecurePointList(pe->GetQEGeom()->GetOrigin(),
						pe->GetQEGeom()->GetDestination());
			}
			++ecIt;
		}
	}
}


}

#endif /* WARPQEMESHFILTER_HXX_ */
