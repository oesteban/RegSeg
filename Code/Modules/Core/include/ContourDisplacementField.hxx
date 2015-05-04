// --------------------------------------------------------------------------
// File:             ContourDisplacementField.hxx
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
// This file is part of RegSeg
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

#ifndef CONTOURDISPLACEMENTFIELD_HXX_
#define CONTOURDISPLACEMENTFIELD_HXX_

#include "ContourDisplacementField.h"

namespace rstk {

template< typename TCoordPrecision, unsigned int VDimension>
ContourDisplacementField< TCoordPrecision, VDimension>::ContourDisplacementField() {
	//this->m_ReferenceMesh = InternalMeshType::New();

}

template< typename TCoordPrecision, unsigned int VDimension>
void
ContourDisplacementField< TCoordPrecision, VDimension>::
PrintSelf( std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf( os, indent );

	os << "* Internal Mesh: " << std::endl;
	//this->m_ReferenceMesh->Print( os, indent.GetNextIndent() );
}

template< typename TCoordPrecision, unsigned int VDimension>
void
ContourDisplacementField< TCoordPrecision, VDimension>::
Initialize() {
	typename Superclass::PointsContainerIterator d_it = this->GetPoints()->Begin();
	VectorType zeroVect = itk::NumericTraits<VectorType>::Zero;
	while (d_it != this->GetPoints()->End()) {
		this->SetPointData( d_it.Index(), zeroVect );
		++d_it;
	}
}

/*
template< typename TCoordPrecision, unsigned int VDimension>
void
ContourDisplacementField< TCoordPrecision, VDimension>::
SetReferenceMesh(InternalMeshType* mesh) {
	typename InternalMeshType::PointsContainerPointer points = mesh->GetPoints();
	typename InternalMeshType::PointsContainerIterator p_it =  points->Begin();

	typename InternalMeshType::PointDataContainerPointer container = mesh->GetPointData();
	typename InternalMeshType::PointDataContainerIterator d_it =     mesh->Begin();

	VectorType zeroVect();
	zeroVect.Fill( 0.0 );
	while (p_it != points->End()) {
		typename Superclass::PointIdentifier id = this->AddPoint( *p_it );
		this->SetPointData( id, zeroVect );
		++p_it;
	}

	typename InternalMeshType::CellsContainerPointer cells = mesh->GetCells();
	typename InternalMeshType::CellsContainerIterator c_it =  cells->Begin();

	while( c_it!= cells->End() ) {
		this->AddFace( )
		++c_it;
	}


}*/

}


#endif /* CONTOURDISPLACEMENTFIELD_HXX_ */
