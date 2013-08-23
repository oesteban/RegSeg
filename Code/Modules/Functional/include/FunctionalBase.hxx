// --------------------------------------------------------------------------
// File:             FunctionalBase.hxx
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

#include "FunctionalBase.h"

#include <iostream>
#include <iomanip>
#include "DisplacementFieldFileWriter.h"

#include <itkMeshFileWriter.h>
#include <itkImageAlgorithm.h>

namespace rstk {


template< typename TReferenceImageType, typename TCoordRepType >
FunctionalBase<TReferenceImageType, TCoordRepType>
::FunctionalBase() {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_FieldInterpolator = FieldInterpolatorType::New();
	this->m_Derivative = FieldType::New();
	this->m_CurrentDisplacementField = FieldType::New();
	this->m_EnergyResampler = DisplacementResamplerType::New();
	this->m_Modified = false;
	this->m_RegionsModified = false;
	this->m_NumberOfContours = 0;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::Initialize() {
	// Set number of control points in the sparse-dense interpolator
	this->m_NumberOfPoints = 0;
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		this->m_NumberOfPoints+= this->m_CurrentContourPosition[contid]->GetNumberOfPoints();
	}
	this->m_FieldInterpolator->SetN(this->m_NumberOfPoints);


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
		PointsVector points;

		// Compute mesh of normals
		NormalFilterPointer normalsFilter = NormalFilterType::New();
		normalsFilter->SetInput( this->m_CurrentContourPosition[contid] );
		normalsFilter->Update();

		ContourPointer normals = normalsFilter->GetOutput();
		outerVect.resize( normals->GetNumberOfPoints() );
		points.resize( normals->GetNumberOfPoints() );

		typename ContourType::PointsContainerConstIterator c_it  = normals->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

		PointType p;
		VectorType v;
		VectorType ni;

		long unsigned int pid;
		while( c_it != c_end ) {
			pid = c_it.Index();
			p = normals->GetPoint(pid);
			points[pid] = p;
			this->m_FieldInterpolator->SetPoint( cpid, p );

			if( contid > 0 ) {
				normals->GetPointData( pid, &ni );
				outerVect[pid] =  interp->Evaluate( p - ni*0.5 );
			}
			++c_it;
			cpid++;
		}
		this->m_ShapePrior.push_back( points );

		if( contid > 0 )
			this->m_OuterList.push_back( outerVect );
	}

	if( this->m_NumberOfContours == 1 ) {
		ContourOuterRegions outerVect;
		outerVect.resize(this->m_CurrentContourPosition[0]->GetNumberOfPoints());
		std::fill( outerVect.begin(), outerVect.end(), 1 );
		this->m_OuterList.push_back( outerVect );
	}

	// Set up grid points in the sparse-dense interpolator
	if( this->m_Derivative.IsNotNull() ) {
		PointType p;
		this->m_NumberOfNodes = this->m_Derivative->GetLargestPossibleRegion().GetNumberOfPixels();
		this->m_FieldInterpolator->SetNumberOfParameters( this->m_NumberOfNodes );

		for ( size_t gid = 0; gid < this->m_NumberOfNodes; gid++ ) {
			this->m_Derivative->TransformIndexToPhysicalPoint( this->m_Derivative->ComputeIndex( gid ), p );
			this->m_FieldInterpolator->SetNode( gid, p );
		}
	} else {
		itkWarningMacro( << "No parametrization (deformation field grid) was defined.");
	}

	typename FieldType::SpacingType sigma = this->m_Derivative->GetSpacing();
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
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContourType* prior ) {
	this->m_CurrentContourPosition.push_back( prior );
	this->m_NumberOfContours++;
    this->m_ROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentMaps.resize( this->m_NumberOfContours+1 );
	return this->m_NumberOfContours-1;
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	this->m_Value = 0.0;

	double normalizer = 1.0;

	for(size_t i = 0; i<Dimension; i++)
		normalizer *= m_ReferenceImage->GetSpacing()[i];

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
	return normalizer*this->m_Value;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::UpdateContour(const typename FunctionalBase<TReferenceImageType, TCoordRepType>::FieldType* newField ) {
	// Copy newField values to interpolator
	const VectorType* NodesBuffer = newField->GetBufferPointer();
	VectorType v;
	for( size_t gpid = 0; gpid < this->m_NumberOfNodes; gpid++ ) {
		v =  *( NodesBuffer + gpid );
		if ( v.GetNorm()>0 ) {
			this->m_FieldInterpolator->SetCoefficient( gpid, v );
		}
	}

	this->m_FieldInterpolator->Interpolate();
	this->m_FieldInterpolator->ComputeNodesData();

	VectorType* fBuffer = this->m_CurrentDisplacementField->GetBufferPointer();

	for( size_t i = 0; i<this->m_NumberOfNodes; i++) {
		*(fBuffer+i) = this->m_FieldInterpolator->GetNodeData(i);
	}

	MeasureType norm;
	MeasureType meanNorm = 0.0;
	VectorType meanDesp;
	meanDesp.Fill(0.0);

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
			desp = this->m_FieldInterpolator->GetPointData( gpid );
			norm = desp.GetNorm();
			// Add vector to the point
			if( norm > 1.0e-3 ) {
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
}

template< typename TReferenceImageType, typename TCoordRepType >
inline bool
FunctionalBase<TReferenceImageType, TCoordRepType>
::IsInside( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType p, typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContinuousIndex& idx) const {
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
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeDerivative() {
	size_t cpid = 0;

	NormalFilterPointer normalsFilter;

#ifndef NDEBUG
	double maxGradient = 0.0;
	double sumGradient = 0.0;
#endif

	VectorType zerov; zerov.Fill(0.0);
	this->m_Derivative->FillBuffer( zerov );

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
			this->m_FieldInterpolator->SetPointData(cpid, ni);
			++c_it;
			cpid++;

#ifndef NDEBUG
			if ( gradient > maxGradient ) maxGradient = gradient;
			sumGradient+=gradient;
#endif
		}
	}

#ifndef NDEBUG
	std::cout << "Gradient: avg=" << ( sumGradient/ cpid ) << ", max=" << maxGradient << "." << std::endl;
	sumGradient = 0.0;
#endif

	// Interpolate sparse velocity field to targetDeformation
	this->m_FieldInterpolator->ComputeWeights();

	VectorType* gmBuffer = this->m_Derivative->GetBufferPointer();
	VectorType v;

	for( size_t gpid = 0; gpid<this->m_NumberOfNodes; gpid++ ) {
		v = this->m_FieldInterpolator->GetNodeWeight(gpid);
		if ( v.GetNorm() > 1.0e-4 ) {
			*( gmBuffer + gpid ) = v; // * (1.0/this->m_NumberOfPoints);
			//this->m_Derivative->SetPixel( this->m_Derivative->ComputeIndex(gpid), v);
#ifndef NDEBUG
			sumGradient+=(*( gmBuffer + gpid )).GetNorm();
			if (v.GetNorm()>maxGradient) maxGradient=(*( gmBuffer + gpid )).GetNorm();
#endif
		}
	}

#ifndef NDEBUG
	std::cout << "Gradient_projected: avg="<< ( sumGradient/ this->m_NumberOfNodes ) << "; max=" << maxGradient << "." << std::endl;
#endif

}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::ROIConstPointer
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetCurrentRegion( size_t idx ) {
	if(this->m_RegionsModified )
		this->ComputeCurrentRegions();

	return this->m_CurrentROIs[idx];
}

template< typename TReferenceImageType, typename TCoordRepType >
const typename FunctionalBase<TReferenceImageType, TCoordRepType>::ProbabilityMapType*
FunctionalBase<TReferenceImageType, TCoordRepType>
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
typename FunctionalBase<TReferenceImageType, TCoordRepType>::ROIConstPointer
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetCurrentRegions() {
	return this->m_CurrentRegions;
}


template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
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
FunctionalBase<TReferenceImageType, TCoordRepType>
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
FunctionalBase<TReferenceImageType, TCoordRepType>
::CopyInformation( const FieldType* field) {
	this->m_Derivative->SetDirection( field->GetDirection() );
	this->m_Derivative->SetOrigin   ( field->GetOrigin() );
	this->m_Derivative->SetSpacing  ( field->GetSpacing() );
	this->m_Derivative->SetRegions  ( field->GetRequestedRegion().GetSize());
	this->m_Derivative->Allocate();

	this->m_CurrentDisplacementField->SetDirection( field->GetDirection() );
	this->m_CurrentDisplacementField->SetOrigin   ( field->GetOrigin() );
	this->m_CurrentDisplacementField->SetSpacing  ( field->GetSpacing() );
	this->m_CurrentDisplacementField->SetRegions  ( field->GetRequestedRegion().GetSize());
	this->m_CurrentDisplacementField->Allocate();
}

}

#endif /* LEVELSETSBASE_HXX_ */
