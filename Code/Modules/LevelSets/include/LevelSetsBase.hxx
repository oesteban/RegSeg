// --------------------------------------------------------------------------
// File:             LevelSetsBase.hxx
// Date:             06/11/2012
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
// This file is part of ACWEReg
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

#ifndef LEVELSETSBASE_HXX_
#define LEVELSETSBASE_HXX_

namespace rstk {



template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::SetShapePrior( typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ContourDeformationType* prior ) {
	this->m_CurrentContourPosition = prior;

	this->m_ContourDeformation = ContourDeformationType::New();

	typename ContourDeformationType::PointsContainerPointer points = prior->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it = points->Begin();
	VectorType zero = itk::NumericTraits<VectorType>::Zero;

	while( u_it != points->End() ) {
		this->m_ContourDeformation->SetPointData(
				this->m_ContourDeformation->AddPoint( u_it.Value() ),zero);
		++u_it;
	}

	typename ContourDeformationType::CellsContainerPointer cells = prior->GetCells();
	typename ContourDeformationType::CellsContainerConstIterator c_it = cells->Begin();

	size_t i = 0;
	while( c_it!=cells->End() ) {
		typename ContourDeformationType::CellType::CellAutoPointer cellCopy;
		c_it.Value()->MakeCopy( cellCopy );
		this->m_ContourDeformation->SetCell( i++ ,cellCopy );
		++c_it;
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::UpdateDeformationField(const typename LevelSetsBase<TReferenceImageType, TCoordRepType>::DeformationFieldType* newField ) {
	// Set-up a linear interpolator for the vector field
	VectorInterpolatorPointer interp = VectorInterpolatorType::New();
	interp->SetInputImage( newField );

	typename ContourDeformationType::PointsContainerPointer points = this->m_ContourDeformation->GetPoints();
	typename ContourDeformationType::PointsContainerIterator p_end = points->End();
	// For all the points in the mesh
	for(typename ContourDeformationType::PointsContainerIterator p_it = points->Begin(); p_it != p_end; ++p_it) {
		PointType currentPoint = p_it.Value();
		// Interpolate the value of the field in the point
		VectorType desp = interp->Evaluate( currentPoint );
		// Add vector to the point
		this->m_CurrentContourPosition->SetPoint( p_it.Index(), currentPoint+desp );
	}
	// TODO this->m_DeformationField = newField;
}

}


#endif /* LEVELSETSBASE_HXX_ */
