// --------------------------------------------------------------------------------------
// File:          SparseDeformationField.txx
// Date:          Apr 4, 2012
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       0.1
// License:       BSD
// Short Summary: 
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Biomedical Image Technology, UPM (BIT-UPM)
// and Signal Processing Lab 5, EPFL (LTS5-EPFL)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//   * Neither the names of the copyright holders, nor the names of
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
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

#ifndef SPARSEDEFORMATIONFIELD_TXX_
#define SPARSEDEFORMATIONFIELD_TXX_

#include "SparseDeformationField.h"

namespace rstk {

template <class TReferenceSurface, class TCoordinatesPrecision >
SparseDeformationField<TReferenceSurface,TCoordinatesPrecision>::
SparseDeformationField() {
	m_Field = SparseVectorFieldType::New();
	m_TransformedPointSet = ReferenceSurfaceType::New();
	m_NumberOfParameters = 0;
}

template <class TReferenceSurface, class TCoordinatesPrecision >
void
SparseDeformationField<TReferenceSurface,TCoordinatesPrecision>::
SetReferenceSurface( const ReferenceSurfaceType * ref ) {
	m_ReferenceSurface = ref;
	m_NumberOfParameters = m_ReferenceSurface->GetNumberOfPoints();
	ParameterType zero;
	zero.Fill( 0.0 );

	m_TransformedPointSet->CopyInformation( m_ReferenceSurface );

	for ( size_t i = 0; i < m_NumberOfParameters; i++) {
		typename ReferenceSurfaceType::PointType point = m_ReferenceSurface->GetPoint(i);
		m_Field->SetPoint( i, point );
		m_TransformedPointSet->SetPoint( i, point );
		this->SetK( i, zero );
	}

	const CellsContainer* cells = m_ReferenceSurface->GetCells( );
	CellsContainerConstIterator itCells = cells->Begin();
	CellsContainerConstIterator itCellsEnd = cells->End();

	size_t i = 0;
    while ( itCells != itCellsEnd )    {
    	typename CellType::CellAutoPointer cellCopy;
        itCells.Value()->MakeCopy( cellCopy );
        m_TransformedPointSet->SetCell( i++, cellCopy );
        ++itCells;
    }
}

template <class TReferenceSurface, class TCoordinatesPrecision >
const typename SparseDeformationField<TReferenceSurface,TCoordinatesPrecision>::ReferenceSurfaceType*
SparseDeformationField<TReferenceSurface,TCoordinatesPrecision>::
GetTransformedSurface() {
	for ( size_t i = 0; i < m_NumberOfParameters; i++) {
		m_TransformedPointSet->SetPoint( i, m_ReferenceSurface->GetPoint(i) + this->GetK(i) );
	}

	return m_TransformedPointSet;
}

}

#endif /* SPARSEDEFORMATIONFIELD_TXX_ */
