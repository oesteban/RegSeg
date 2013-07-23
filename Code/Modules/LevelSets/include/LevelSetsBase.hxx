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
#include "DisplacementFieldFileWriter.h"

#include <itkMeshFileWriter.h>
#include <itkImageAlgorithm.h>

namespace rstk {


template< typename TReferenceImageType, typename TCoordRepType >
LevelSetsBase<TReferenceImageType, TCoordRepType>
::LevelSetsBase() {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_FieldInterpolator = FieldInterpolatorType::New();
	this->m_GradientMap = FieldType::New();
	this->m_EnergyResampler = DisplacementResamplerType::New();
	this->m_Modified = false;
	this->m_RegionsModified = false;
	this->m_NumberOfContours = 0;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::Initialize() {
	// Set number of control points in the sparse-dense interpolator
	size_t nControlPoints = 0;
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		nControlPoints+= this->m_CurrentContourPosition[contid]->GetNumberOfPoints();
	}
	this->m_FieldInterpolator->SetN(nControlPoints);


	// Initialize corresponding ROI /////////////////////////////
	// Check that high-res reference sampling grid has been initialized
	if ( this->m_ReferenceSamplingGrid.IsNull() ) {
			this->InitializeSamplingGrid();
			this->m_CurrentRegions = ROIType::New();
			this->m_CurrentRegions->SetSpacing(   this->m_ReferenceSamplingGrid->GetSpacing() );
			this->m_CurrentRegions->SetDirection( this->m_ReferenceSamplingGrid->GetDirection() );
			this->m_CurrentRegions->SetOrigin(    this->m_ReferenceSamplingGrid->GetOrigin() );
			this->m_CurrentRegions->SetRegions(   this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
			this->m_CurrentRegions->Allocate();
	}

	// Compute and set regions
	this->ComputeCurrentRegions();
	for( size_t id = 0; id < m_ROIs.size(); id++) {
		this->m_ROIs[id] = this->m_CurrentROIs[id];
	}

	// Set up ROI interpolator
	typename ROIInterpolatorType::Pointer interp = ROIInterpolatorType::New();
	interp->SetInputImage( this->m_CurrentRegions );


	// Set up outer regions AND control points in the interpolator
	size_t cpid = 0;
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		ContourOuterRegions outerVect;
		ControlPointsVector controlPoints;

		// Compute mesh of normals
		NormalFilterPointer normalsFilter = NormalFilterType::New();
		normalsFilter->SetInput( this->m_CurrentContourPosition[contid] );
		normalsFilter->Update();

		ContourPointer normals = normalsFilter->GetOutput();
		outerVect.resize( normals->GetNumberOfPoints() );
		controlPoints.resize( normals->GetNumberOfPoints() );

		typename ContourType::PointsContainerConstIterator c_it  = normals->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

		PointType p;
		VectorType v;
		VectorType ni;

		long unsigned int pid;
		while( c_it != c_end ) {
			pid = c_it.Index();
			p = normals->GetPoint(pid);
			controlPoints[pid] = p;
			this->m_FieldInterpolator->SetControlPoint( cpid, p );

			if( this->m_NumberOfContours>1 )
				normals->GetPointData( pid, &ni );
				outerVect[pid] =  interp->Evaluate( p - ni*0.5 );
			++c_it;
			cpid++;
		}
		this->m_ShapePrior.push_back( controlPoints );

		if( this->m_NumberOfContours>1 )
			this->m_OuterList.push_back( outerVect );
	}

	if( this->m_NumberOfContours == 1 ) {
		ContourOuterRegions outerVect;
		outerVect.resize(this->m_CurrentContourPosition[0]->GetNumberOfPoints());
		std::fill( outerVect.begin(), outerVect.end(), 1 );
		this->m_OuterList.push_back( outerVect );
	}



	// Set up grid points in the sparse-dense interpolator
	if( this->m_GradientMap.IsNotNull() ) {
		PointType p;
		size_t nGridPoints = this->m_GradientMap->GetLargestPossibleRegion().GetNumberOfPixels();
		this->m_FieldInterpolator->SetK( nGridPoints );

		for ( size_t gid = 0; gid < nGridPoints; gid++ ) {
			this->m_GradientMap->TransformIndexToPhysicalPoint( this->m_GradientMap->ComputeIndex( gid ), p );
			this->m_FieldInterpolator->SetGridPoint( gid, p );
		}
	} else {
		itkWarningMacro( << "No parametrization (deformation field grid) was defined.");
	}

	typename FieldType::SpacingType sigma = this->m_GradientMap->GetSpacing();
	this->m_FieldInterpolator->SetSigma( sigma*2.0 );

#ifndef DNDEBUG
	for( size_t id = 0; id < m_ROIs.size(); id++) {
		typedef itk::ImageFileWriter< ROIType > ROIWriter;
		typename ROIWriter::Pointer w = ROIWriter::New();
		w->SetInput( this->m_ROIs[id] );
		std::stringstream ss;
		ss << "roi_" << std::setfill( '0' ) << std::setw(2) << id << ".nii.gz";
		w->SetFileName( ss.str().c_str() );
		w->Update();
	}
	typedef itk::ImageFileWriter< ROIType > ROIWriter;
	typename ROIWriter::Pointer w = ROIWriter::New();
	w->SetInput( this->m_CurrentRegions );
	w->SetFileName( "regions.nii.gz" );
	w->Update();
#endif
}


template< typename TReferenceImageType, typename TCoordRepType >
size_t
LevelSetsBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ContourType* prior ) {
	this->m_CurrentContourPosition.push_back( prior );
	this->m_NumberOfContours++;
    this->m_ROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentMaps.resize( this->m_NumberOfContours+1 );
	return this->m_NumberOfContours-1;
}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::MeasureType
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	this->m_Value = 0.0;

	for( size_t roi = 0; roi < m_ROIs.size(); roi++ ) {
		ProbabilityMapConstPointer roipm = this->GetCurrentMap( roi );

		const typename ProbabilityMapType::PixelType* roiBuffer = roipm->GetBufferPointer();

		size_t nPix = roipm->GetLargestPossibleRegion().GetNumberOfPixels();
		ReferencePointType pos, targetPos;
		typename ProbabilityMapType::PixelType w;

		for( size_t i = 0; i < nPix; i++) {
			w = *( roiBuffer + i );
			if ( w > 0.0 ) {
				roipm->TransformIndexToPhysicalPoint( roipm->ComputeIndex(i), pos);
				this->m_Value +=  w * this->GetEnergyAtPoint( targetPos, roi );
			}
		}
	}
	return this->m_Value;
}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::MeasureType
LevelSetsBase<TReferenceImageType, TCoordRepType>
::UpdateContour(const typename LevelSetsBase<TReferenceImageType, TCoordRepType>::FieldType* newField ) {
	// Copy newField values to interpolator
	size_t nGridPoints = newField->GetLargestPossibleRegion().GetNumberOfPixels();
	const VectorType* gridPointsBuffer = newField->GetBufferPointer();
	VectorType v;
	MeasureType maxNorm = 0.0;
	for( size_t gpid = 0; gpid < nGridPoints; gpid++ ) {
		v =  *( gridPointsBuffer + gpid );
		if ( v.GetNorm()>0 ) {
			this->m_FieldInterpolator->SetGridPointData( gpid, v );

			if( v.GetNorm()>maxNorm) maxNorm=v.GetNorm();
		}
	}

#ifndef NDEBUG
	std::cout << "Field maxNorm=" << maxNorm << "mm." << std::endl;
#endif


	this->m_FieldInterpolator->ComputeControlPoints();

	MeasureType norm;
	MeasureType meanNorm = 0.0;
	maxNorm = 0.0;
	VectorType meanDesp;
	meanDesp.Fill(0.0);

	// Set-up a linear interpolator for the vector field
	VectorInterpolatorPointer interp = VectorInterpolatorType::New();
	interp->SetInputImage( newField );

	ContinuousIndex point_idx;
	size_t changed = 0;

	size_t gpid = 0;
	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++ ) {
		typename ContourType::PointsContainerPointer curr_points = this->m_CurrentContourPosition[contid]->GetPoints();
		typename ContourType::PointsContainerIterator p_it = curr_points->Begin();
		typename ContourType::PointsContainerIterator p_end = curr_points->End();

		// For all the points in the mesh
		PointType p;
		VectorType desp, desp_back;
		typename ContourType::PointType currentPoint,newPoint;
		long unsigned int pid;

		while ( p_it != p_end ) {
			currentPoint = p_it.Value();
			pid = p_it.Index();

			// Interpolate the value of the field in the point
			desp_back = interp->Evaluate( currentPoint );
			desp = this->m_FieldInterpolator->GetGridPointData( gpid );
			norm = desp.GetNorm();
			// Add vector to the point
			if( norm > 1.0e-3 ) {
				if ( norm > maxNorm ) maxNorm = norm;
				meanDesp += desp;
				meanNorm += norm;
				p = this->m_ShapePrior[contid][pid] + desp;
				newPoint.SetPoint( p );
				newPoint.SetEdge( currentPoint.GetEdge() );

				if( ! this->IsInside(p,point_idx) ) {
					itkExceptionMacro( << "Contour is outside image regions after update.\n" <<
						"\tMoving vertex [" << contid << "," << pid << "] from " << this->m_ShapePrior[contid][pid] << " to " << p << " norm="  << desp.GetNorm() << "mm.\n");
				}

				this->m_CurrentContourPosition[contid]->SetPoint( p_it.Index(), newPoint );
				changed++;
				gpid++;
			}

			++p_it;
		}
	}

	this->m_RegionsModified = (changed>0);

#ifndef NDEBUG
	std::cout << "MeanNorm=" << (meanNorm/changed) << "mm.; maxNorm=" << maxNorm << "mm.;" << " meanDesp=" << (meanDesp/changed) << std::endl;
#endif

	return (meanNorm/changed);
}

template< typename TReferenceImageType, typename TCoordRepType >
inline bool
LevelSetsBase<TReferenceImageType, TCoordRepType>
::IsInside( const typename LevelSetsBase<TReferenceImageType, TCoordRepType>::PointType p, typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ContinuousIndex& idx) const {
	bool isInside = this->m_ReferenceImage->TransformPhysicalPointToContinuousIndex( p, idx );
#ifndef NDEBUG
	if(!isInside) {
		typename FieldType::PointType origin, end;
		typename FieldType::IndexType tmp_idx;
		typename FieldType::SizeType size = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
		tmp_idx.Fill(0);
		this->m_ReferenceImage->TransformIndexToPhysicalPoint( tmp_idx, origin );
		for ( size_t dim = 0; dim<FieldType::ImageDimension; dim++)  tmp_idx[dim]= size[dim]-1;
		this->m_ReferenceImage->TransformIndexToPhysicalPoint( tmp_idx, end);
		itkWarningMacro( << "Point p=[" << p << "] is outside image extents (" << origin << ", " << end << ").");
	}
#endif
	return isInside;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::ComputeGradient() {
	size_t cpid = 0;

	NormalFilterPointer normalsFilter;

#ifndef NDEBUG
	double maxGradient = 0.0;
	double sumGradient = 0.0;
#endif

	VectorType zerov; zerov.Fill(0.0);
	this->m_GradientMap->FillBuffer( zerov );

	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++) {
		// Compute mesh of normals
		normalsFilter = NormalFilterType::New();
		normalsFilter->SetInput( this->m_CurrentContourPosition[contid] );
		normalsFilter->Update();
		ContourPointer normals = normalsFilter->GetOutput();

#ifndef NDEBUG
		typedef itk::MeshFileWriter< ContourType >     MeshWriterType;
		typename MeshWriterType::Pointer w = MeshWriterType::New();
		w->SetInput( this->m_CurrentContourPosition[contid] );
		std::stringstream ss;
		ss << "normals_" << std::setfill('0') << std::setw(2) << contid << ".vtk";
		w->SetFileName( ss.str() );
		w->Update();
#endif

		typename ContourType::PointsContainerConstIterator c_it = normals->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

		PointValueType gradient;
		PointType  ci, ci_prime;
		VectorType ni;
		typename ContourType::PointIdentifier pid;
		ReferencePixelType dist[2];
		size_t outer_contid;

		// for all node in mesh
		while (c_it!=c_end) {
			pid = c_it.Index();
			outer_contid = this->m_OuterList[contid][pid];

			if ( contid != outer_contid ) {
				ci_prime = c_it.Value();
				gradient =  this->GetEnergyAtPoint( ci_prime, outer_contid ) - this->GetEnergyAtPoint( ci_prime, contid );
				assert( !std::isnan(gradient) );
			} else {
				gradient = 0.0;
			}

			// project to normal, updating transform
			normals->GetPointData( pid, &ni );         // Normal ni in point c'_i
			ni*= gradient;
			normals->SetPointData( pid, ni );
			this->m_FieldInterpolator->SetControlPointData(cpid, ni);
			++c_it;
			cpid++;

#ifndef NDEBUG
			if ( gradient > maxGradient ) maxGradient = gradient;
			sumGradient+=gradient;
#endif
		}
	}

#ifndef NDEBUG
	std::cout << "max_gradient=" << maxGradient << ", avg=" << ( sumGradient/ cpid ) << "." << std::endl;
#endif

	// Interpolate sparse velocity field to targetDeformation
	this->m_FieldInterpolator->ComputeGridPoints();

	VectorType* gmBuffer = this->m_GradientMap->GetBufferPointer();
	size_t nPoints = this->m_GradientMap->GetLargestPossibleRegion().GetNumberOfPixels();
	VectorType v;

	sumGradient = 0.0;
	for( size_t gpid = 0; gpid<nPoints; gpid++ ) {
		v = this->m_FieldInterpolator->GetGridPointData( gpid );
		if ( v.GetNorm() > 0 ) {
			this->m_GradientMap->SetPixel( this->m_GradientMap->ComputeIndex(gpid), v);
			sumGradient+=v.GetNorm();
			if (v.GetNorm()>maxGradient) maxGradient=v.GetNorm();
		}
	}

#ifndef NDEBUG
	std::cout << "max_interpolated_gradient=" << maxGradient << ", avg=" << ( sumGradient/ nPoints ) << "." << std::endl;
#endif


}

template< typename TReferenceImageType, typename TCoordRepType >
typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ROIConstPointer
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetCurrentRegion( size_t idx ) {
	if(this->m_RegionsModified )
		this->ComputeCurrentRegions();

	return this->m_CurrentROIs[idx];
}

template< typename TReferenceImageType, typename TCoordRepType >
const typename LevelSetsBase<TReferenceImageType, TCoordRepType>::ProbabilityMapType*
LevelSetsBase<TReferenceImageType, TCoordRepType>
::GetCurrentMap( size_t idx ) {
	if(this->m_RegionsModified ) {
		this->ComputeCurrentRegions();
	}

	if( this->m_CurrentMaps[idx].IsNull() ) {
		this->m_CurrentMaps[idx] = ProbabilityMapType::New();
		this->m_CurrentMaps[idx]->SetRegions(      this->m_ReferenceImage->GetLargestPossibleRegion().GetSize() );
		this->m_CurrentMaps[idx]->SetOrigin(    this->m_ReferenceImage->GetOrigin() );
		this->m_CurrentMaps[idx]->SetDirection( this->m_ReferenceImage->GetDirection() );
		this->m_CurrentMaps[idx]->SetSpacing(   this->m_ReferenceImage->GetSpacing() );
		this->m_CurrentMaps[idx]->Allocate();
		this->m_CurrentMaps[idx]->FillBuffer( 0.0 );
	}

	// Resample to reference image resolution
	ResampleROIFilterPointer resampleFilter = ResampleROIFilterType::New();
	resampleFilter->SetInput( this->m_CurrentROIs[idx] );
	resampleFilter->SetSize( this->m_ReferenceImage->GetLargestPossibleRegion().GetSize() );
	resampleFilter->SetOutputOrigin(    this->m_ReferenceImage->GetOrigin() );
	resampleFilter->SetOutputSpacing(   this->m_ReferenceImage->GetSpacing() );
	resampleFilter->SetOutputDirection( this->m_ReferenceImage->GetDirection() );
	resampleFilter->SetDefaultPixelValue( 0.0 );
	resampleFilter->Update();
	ProbabilityMapPointer tpm = resampleFilter->GetOutput();

	itk::ImageAlgorithm::Copy<ProbabilityMapType,ProbabilityMapType>(
			tpm, this->m_CurrentMaps[idx],
			tpm->GetLargestPossibleRegion(),
			this->m_CurrentMaps[idx]->GetLargestPossibleRegion()
	);

	return this->m_CurrentMaps[idx].GetPointer();
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

	this->m_RegionsModified = false;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::InitializeSamplingGrid() {
	this->m_ReferenceSamplingGrid = FieldType::New();
	typename ReferenceImageType::SpacingType sp = this->m_ReferenceImage->GetSpacing();
	double spacing = itk::NumericTraits< double >::max();
	double factor = 0.25;

	for (size_t i = 0; i<Dimension; i++ ){
		if (sp[i] < spacing )
			spacing = sp[i];
	}

	sp.Fill( spacing * factor );

	PointType origin = this->m_ReferenceImage->GetOrigin();
	typename ReferenceImageType::SizeType size = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
	typename itk::ContinuousIndex<double, Dimension> idx;
	for (size_t i = 0; i<Dimension; i++ ){
		idx[i]=size[i];
	}

	PointType end;
	this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( idx, end );

	for (size_t i = 0; i<Dimension; i++ ){
		size[i]= (unsigned int) ( fabs( (end[i]-origin[i])/sp[i] ) );
	}

	this->m_ReferenceSamplingGrid->SetOrigin( origin );
	this->m_ReferenceSamplingGrid->SetDirection( this->m_ReferenceImage->GetDirection() );
	this->m_ReferenceSamplingGrid->SetRegions( size );
	this->m_ReferenceSamplingGrid->SetSpacing( sp );
	this->m_ReferenceSamplingGrid->Allocate();

#ifndef NDEBUG
	typedef itk::VectorResampleImageFilter< ReferenceImageType, ReferenceImageType > Int;
	typename Int::Pointer intp = Int::New();
	intp->SetInput( this->m_ReferenceImage );
	intp->SetOutputSpacing( this->m_ReferenceSamplingGrid->GetSpacing() );
	intp->SetOutputDirection( this->m_ReferenceSamplingGrid->GetDirection() );
	intp->SetOutputOrigin( this->m_ReferenceSamplingGrid->GetOrigin() );
	intp->SetSize( this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
	intp->Update();

	typedef itk::ImageFileWriter< ReferenceImageType > W;
	typename W::Pointer w = W::New();
	w->SetInput( intp->GetOutput() );
	w->SetFileName( "ReferenceSamplingGridTest.nii.gz" );
	w->Update();
#endif
}

template< typename TReferenceImageType, typename TCoordRepType >
void
LevelSetsBase<TReferenceImageType, TCoordRepType>
::CopyInformation( const FieldType* field) {
	this->m_GradientMap->SetDirection( field->GetDirection() );
	this->m_GradientMap->SetOrigin   ( field->GetOrigin() );
	this->m_GradientMap->SetSpacing  ( field->GetSpacing() );
	this->m_GradientMap->SetRegions  ( field->GetRequestedRegion().GetSize());
	this->m_GradientMap->Allocate();
}

}

#endif /* LEVELSETSBASE_HXX_ */
