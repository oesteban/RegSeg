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

#include <iostream>
#include <iomanip>

namespace rstk {


template< typename TReferenceImageType, typename TCoordRepType >
LevelSetsBase<TReferenceImageType, TCoordRepType>
::LevelSetsBase() {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_SparseToDenseResampler = SparseToDenseFieldResampleType::New();
	this->m_EnergyResampler = DisplacementResamplerType::New();
	this->m_Modified = false;
	this->m_RegionsModified = false;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::Initialize() {
	// Initialize corresponding ROI /////////////////////////////
	// 1. Check that high-res reference sampling grid has been initialized
	if ( this->m_ReferenceSamplingGrid.IsNull() ) {
			this->InitializeSamplingGrid();
			this->m_CurrentRegions = ROIType::New();
			this->m_CurrentRegions->SetSpacing(   this->m_ReferenceSamplingGrid->GetSpacing() );
			this->m_CurrentRegions->SetDirection( this->m_ReferenceSamplingGrid->GetDirection() );
			this->m_CurrentRegions->SetOrigin(    this->m_ReferenceSamplingGrid->GetOrigin() );
			this->m_CurrentRegions->SetRegions(      this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
			this->m_CurrentRegions->Allocate();
	}

	this->ComputeCurrentRegions();

	for( size_t id = 0; id < m_ROIs.size(); id++) {
		this->m_ROIs[id] = this->m_CurrentROIs[id];
#ifndef DNDEBUG
		typedef itk::ImageFileWriter< ROIType > ROIWriter;
		typename ROIWriter::Pointer w = ROIWriter::New();
		w->SetInput( this->m_ROIs[id] );
		std::stringstream ss;
		ss << "roi_" << std::setfill( '0' ) << std::setw(2) << id << ".nii.gz";
		w->SetFileName( ss.str().c_str() );
		w->Update();
#endif
	}

#ifndef DNDEBUG
	typedef itk::ImageFileWriter< ROIType > ROIWriter;
	typename ROIWriter::Pointer w = ROIWriter::New();
	w->SetInput( this->m_CurrentRegions );
	w->SetFileName( "regions.nii.gz" );
	w->Update();
#endif

	//
	// Assign the outer region to each region
	//

	// 1. Set up ROI interpolator
	typename ROIInterpolatorType::Pointer interp = ROIInterpolatorType::New();
	interp->SetInputImage( this->m_CurrentRegions );

	for ( size_t id = 0; id < this->m_CurrentContourPosition.size(); id ++) {
		ContourOuterRegions outerVect;

		// Compute mesh of normals
		NormalFilterPointer normFilter = NormalFilterType::New();
		normFilter->SetInput( this->m_CurrentContourPosition[id] );
		normFilter->Update();
		ContourDeformationPointer normals = normFilter->GetOutput();

		typename ContourDeformationType::PointsContainerConstIterator c_it     = this->m_ShapePrior[id]->GetPoints()->Begin();
		typename ContourDeformationType::PointsContainerConstIterator c_end    = this->m_ShapePrior[id]->GetPoints()->End();

		PixelPointType p;
		itk::Vector<double, 3u> norm;
		long unsigned int idx;
		while( c_it != c_end ) {
			p = c_it.Value();
			idx = c_it.Index();
			normals->GetPointData( idx, &norm );
			p = p + norm;
			outerVect.push_back(  interp->Evaluate( p ) );
			++c_it;
		}

		this->m_OuterList.push_back( outerVect );
	}

	this->m_Modified = false;
}


template< typename TReferenceImageType, typename TCoordRepType >
size_t
LevelSetsBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ContourDeformationType* prior ) {
	this->m_CurrentContourPosition.push_back( prior );
	/*
	// TODO Check that contour copier now works.

	this->m_ContourCopier = ContourCopyType::New();
	m_ContourCopier->SetInput( prior );
	this->m_ContourCopier->Update();
	this->m_ShapePrior = m_ContourCopier->GetOutput();*/

    this->m_ShapePrior.push_back( ContourDeformationType::New() );
    this->m_ROIs.resize( this->m_ShapePrior.size()+1 );
    this->m_CurrentROIs.resize( this->m_ShapePrior.size()+1 );

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

	return this->m_CurrentContourPosition.size()-1;

	// TODO Check overlap (regions are disjoint see issue #76).

}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::MeasureType
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	this->m_Value = 0.0;

	size_t nPix = this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetNumberOfPixels();

	for( size_t roi = 0; roi < m_ROIs.size(); roi++ ) {
		ROIConstPointer mask = this->GetCurrentRegion(roi);
		const unsigned char * roiBuffer = mask->GetBufferPointer();

		PixelPointType pos, targetPos;
		for( size_t i = 0; i < nPix; i++) {
			if ( *( roiBuffer + i ) > 0.0 ) {
				mask->TransformIndexToPhysicalPoint( mask->ComputeIndex(i), pos);
				this->m_Value +=  this->GetEnergyAtPoint( targetPos, roi );
			}
		}
	}
	return this->m_Value;
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
	size_t changed = 0;

	for( size_t cont = 0; cont < this->m_ShapePrior.size(); cont++ ) {
		typename ContourDeformationType::PointsContainerPointer curr_points = this->m_CurrentContourPosition[cont]->GetPoints();
		typename ContourDeformationType::PointsContainerIterator p_end = curr_points->End();
		typename ContourDeformationType::PointsContainerPointer shape_points = this->m_ShapePrior[cont]->GetPoints();
		typename ContourDeformationType::PointsContainerIterator shape_it = shape_points->Begin();

		// For all the points in the mesh
		PointType currentPoint,newPoint;
		PixelPointType p;
		for(typename ContourDeformationType::PointsContainerIterator p_it = curr_points->Begin(); p_it != p_end; ++p_it, ++shape_it) {
			currentPoint = p_it.Value();
			// Interpolate the value of the field in the point
			VectorType desp = interp->Evaluate( currentPoint );
			// Add vector to the point
			if( desp.GetNorm()>0 ) {
				p = shape_it.Value() + desp;
				newPoint.SetPoint( p );
				newPoint.SetEdge( currentPoint.GetEdge() );

				if(! newField->TransformPhysicalPointToContinuousIndex( p, point_idx ) ) {
					typename DeformationFieldType::PointType origin, end;
					typename DeformationFieldType::IndexType tmp_idx;
					typename DeformationFieldType::SizeType size = newField->GetLargestPossibleRegion().GetSize();
					tmp_idx.Fill(0);
					newField->TransformIndexToPhysicalPoint( tmp_idx, origin );
					for ( size_t dim = 0; dim<DeformationFieldType::ImageDimension; dim++)  tmp_idx[dim]= size[dim]-1;
					newField->TransformIndexToPhysicalPoint( tmp_idx, end);
					itkExceptionMacro( << "Contour is outside image regions after update: " << std::endl <<
							"\tMoving vertex [" << shape_it.Index() << "] from " << shape_it.Value() << " to " << newPoint << " norm="  << desp.GetNorm() << "mm." <<
							"\tImage Extents=(" << origin << ", " << end << ").");
				}

				this->m_CurrentContourPosition[cont]->SetPoint( p_it.Index(), newPoint );
				changed++;
			}
		}
	}

	this->m_Modified = (changed>0);
}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::DeformationFieldPointer
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetLevelSetsMap( LevelSetsBase<TReferenceImageType, TCoordRepType>::DeformationFieldType* levelSetMap) {
#ifndef NDEBUG
	double maxLS = 0.0;
#endif

	for( size_t cont = 0; cont < this->m_CurrentContourPosition.size(); cont++) {
		// Compute mesh of normals
		NormalFilterPointer normFilter = NormalFilterType::New();
		normFilter->SetInput( this->m_CurrentContourPosition[cont] );
		normFilter->Update();
		ContourDeformationPointer normals = normFilter->GetOutput();

		typename ContourDeformationType::PointsContainerConstIterator c_it = normals->GetPoints()->Begin();
		typename ContourDeformationType::PointsContainerConstIterator c_end = normals->GetPoints()->End();


		PointValueType levelSet;
		PointType  ci;          //
		PixelPointType  ci_prime;
		PixelType  fi;          // Feature on ci_prime
		typename ContourDeformationType::PixelType ni;
		typename ContourDeformationType::PointIdentifier idx;
		PixelType dist[2];
		size_t outer_cont;

		// for all node in mesh
		while (c_it!=c_end) {
			idx = c_it.Index();
			outer_cont = this->m_OuterList[cont][idx];
			ci_prime = c_it.Value();

			levelSet = - this->GetEnergyAtPoint( ci_prime, cont )
					   + this->GetEnergyAtPoint( ci_prime, outer_cont );

			assert( !std::isnan(levelSet) );
			// project to normal, updating transform
			normals->GetPointData( idx, &ni );         // Normal ni in point c'_i
			ni*=levelSet;
			normals->SetPointData( idx, ni );
			++c_it;
#ifndef NDEBUG
			if ( fabs( levelSet ) > fabs( maxLS ) )
				maxLS = levelSet;
#endif
		}
		this->m_SparseToDenseResampler->SetInput( cont, normals );
	}

#ifndef NDEBUG
	std::cout << "Shape Gradient maximum update value = " << maxLS << std::endl;
#endif

	// Interpolate sparse velocity field to targetDeformation
	this->m_SparseToDenseResampler->Update();
	return this->m_SparseToDenseResampler->GetOutput();

}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ROIConstPointer
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetCurrentRegion( size_t idx ) {
	if(this->m_Modified )
		this->ComputeCurrentRegions();

	return this->m_CurrentROIs[idx];
}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ROIConstPointer
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetCurrentRegions() {
	return this->m_CurrentRegions;
}


template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::ComputeCurrentRegions() {
	ROIPixelType unassigned = itk::NumericTraits< ROIPixelType >::max();
	this->m_CurrentRegions->FillBuffer(unassigned);
	size_t nPix = this->m_CurrentRegions->GetLargestPossibleRegion().GetNumberOfPixels();


	ROIPixelType* regionsBuffer = this->m_CurrentRegions->GetBufferPointer();

	for (ROIPixelType idx = 0; idx < this->m_CurrentROIs.size(); idx++ ) {
		ROIPointer tempROI;

		if ( idx < this->m_CurrentROIs.size() - 1 ) {
			BinarizeMeshFilterPointer meshFilter = BinarizeMeshFilterType::New();
			meshFilter->SetSpacing(   this->m_ReferenceSamplingGrid->GetSpacing() );
			meshFilter->SetDirection( this->m_ReferenceSamplingGrid->GetDirection() );
			meshFilter->SetOrigin(    this->m_ReferenceSamplingGrid->GetOrigin() );
			meshFilter->SetSize(      this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
			meshFilter->SetInput(     this->m_CurrentContourPosition[idx]);
			meshFilter->Update();
			tempROI = meshFilter->GetOutput();
		} else {
			tempROI = ROIType::New();
			tempROI->SetSpacing(   this->m_ReferenceSamplingGrid->GetSpacing() );
			tempROI->SetDirection( this->m_ReferenceSamplingGrid->GetDirection() );
			tempROI->SetOrigin(    this->m_ReferenceSamplingGrid->GetOrigin() );
			tempROI->SetRegions(      this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
			tempROI->Allocate();
			tempROI->FillBuffer( 1 );
		}

		ROIPixelType* roiBuffer = tempROI->GetBufferPointer();

		for( size_t pix = 0; pix < nPix; pix++ ) {
			if( *(regionsBuffer+pix) == unassigned && *( roiBuffer + pix )==1 ) {
				*(regionsBuffer+pix) = idx;
			} else {
				*( roiBuffer + pix ) = 0;
			}
		}

		this->m_CurrentROIs[idx] = tempROI;
	}
}

}

#endif /* LEVELSETSBASE_HXX_ */
