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

#include "LevelSetsBase.h"

namespace rstk {


template< typename TReferenceImageType, typename TCoordRepType >
LevelSetsBase<TReferenceImageType, TCoordRepType>
::LevelSetsBase() {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_SparseToDenseResampler = SparseToDenseFieldResampleType::New();
}


template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ContourDeformationType* prior ) {
	this->m_CurrentContourPosition.push_back( prior );
	/*
	this->m_ContourCopier = ContourCopyType::New();
	m_ContourCopier->SetInput( prior );
	this->m_ContourCopier->Update();
	this->m_ShapePrior = m_ContourCopier->GetOutput();*/

    this->m_ShapePrior.push_back( ContourDeformationType::New() );
    ContourDeformationPointer shapePrior = this->m_ShapePrior.back();

	typename ContourDeformationType::PointsContainerConstIterator u_it = prior->GetPoints()->Begin();
    typename ContourDeformationType::PointsContainerConstIterator u_end = prior->GetPoints()->End();
	VectorType zero = itk::NumericTraits<VectorType>::Zero;

	PointType p, newP;
	while( u_it != u_end ) {
		p = u_it.Value();
		newP.SetPoint( p );
		newP.SetEdge( p.GetEdge() );
		shapePrior->SetPointData( shapePrior->AddPoint( newP ), zero);
		++u_it;
	}

	typename ContourDeformationType::CellsContainerConstIterator c_it  = prior->GetCells()->Begin();
	typename ContourDeformationType::CellsContainerConstIterator c_end = prior->GetCells()->End();

	size_t i = 0;
	while( c_it!=c_end ) {
		typename ContourDeformationType::CellType::CellAutoPointer cellCopy;
		c_it.Value()->MakeCopy( cellCopy );
		shapePrior->SetCell( i++ ,cellCopy );
		++c_it;
	}

	this->m_SparseToDenseResampler->AddControlPoints( shapePrior );
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::UpdateDeformationField(const typename LevelSetsBase<TReferenceImageType, TCoordRepType>::DeformationFieldType* newField ) {
	/*
	if ( this->m_ContourUpdater.IsNull() ) {
		this->m_ContourUpdater = WarpContourType::New();
		this->m_ContourUpdater->SetInput(this->m_ShapePrior);
		this->m_CurrentContourPosition = this->m_ContourUpdater->GetOutput();
	}

	this->m_ContourUpdater->SetDisplacementField( newField );
	this->m_ContourUpdater->Update();
	*/

	// Set-up a linear interpolator for the vector field
	VectorInterpolatorPointer interp = VectorInterpolatorType::New();
	interp->SetInputImage( newField );

	ContinuousIndex point_idx;
	for( size_t cont = 0; cont < this->m_ShapePrior.size(); cont++ ) {
		typename ContourDeformationType::PointsContainerPointer curr_points = this->m_CurrentContourPosition[cont]->GetPoints();
		typename ContourDeformationType::PointsContainerIterator p_end = curr_points->End();
		typename ContourDeformationType::PointsContainerPointer shape_points = this->m_ShapePrior[cont]->GetPoints();
		typename ContourDeformationType::PointsContainerIterator shape_it = shape_points->Begin();

		// For all the points in the mesh
		PointType currentPoint,newPoint;
		for(typename ContourDeformationType::PointsContainerIterator p_it = curr_points->Begin(); p_it != p_end; ++p_it, ++shape_it) {
			currentPoint = p_it.Value();
			// Interpolate the value of the field in the point
			VectorType desp = interp->Evaluate( currentPoint );
			// Add vector to the point
			if( desp.GetNorm()>0 ) {
				newPoint.SetPoint( shape_it.Value() + desp );
				newPoint.SetEdge( currentPoint.GetEdge() );

				if(! newField->TransformPhysicalPointToContinuousIndex( newPoint, point_idx ) ) {
					typename DeformationFieldType::PointType origin, end;
					typename DeformationFieldType::IndexType tmp_idx;
					typename DeformationFieldType::SizeType size = newField->GetLargestPossibleRegion().GetSize();
					tmp_idx.Fill(0);
					newField->TransformIndexToPhysicalPoint( tmp_idx, origin );
					for ( size_t dim = 0; dim<DeformationFieldType::ImageDimension; dim++)  tmp_idx[dim]= size[dim]-1;
					newField->TransformIndexToPhysicalPoint( tmp_idx, end);
					itkExceptionMacro( << "Contour is outside image regions after update: vertex " << newPoint << " is outside image extents (" << origin << ", " << end << ").");
				}

				this->m_CurrentContourPosition[cont]->SetPoint( p_it.Index(), newPoint );
			}
		}
	}
}

}


#endif /* LEVELSETSBASE_HXX_ */
