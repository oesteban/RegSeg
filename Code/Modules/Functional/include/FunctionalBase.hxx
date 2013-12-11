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
#include <math.h>
#include <numeric>
#include <vnl/vnl_random.h>
#include "DisplacementFieldFileWriter.h"

#include <itkMeshFileWriter.h>
#include <itkImageAlgorithm.h>
#include <itkOrientImageFilter.h>
#include <itkContinuousIndex.h>

namespace rstk {


template< typename TReferenceImageType, typename TCoordRepType >
FunctionalBase<TReferenceImageType, TCoordRepType>
::FunctionalBase() {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_FieldInterpolator = FieldInterpolatorType::New();
	this->m_Derivative = FieldType::New();
	this->m_CurrentDisplacementField = FieldType::New();
	this->m_EnergyResampler = DisplacementResamplerType::New();
	this->m_EnergyUpdated = false;
	this->m_RegionsUpdated = false;
	this->m_NumberOfContours = 0;
	this->m_SamplingFactor = 4;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::Initialize() {
	this->InitializeCurrentContours();

	// Initialize corresponding ROI /////////////////////////////
	// Check that high-res reference sampling grid has been initialized
	if ( this->m_ReferenceSamplingGrid.IsNull() ) {
			this->InitializeSamplingGrid();
	}

	// Compute and set regions in m_ROIs
	this->ComputeCurrentRegions();
	for( size_t id = 0; id < m_ROIs.size(); id++) {
		this->m_ROIs[id] = this->m_CurrentROIs[id];
	}

	// Compute the outer region in each vertex
	this->ComputeOuterRegions();

	// Set up grid points in the sparse-dense interpolator
	this->InitializeInterpolatorGrid();
}


template< typename TReferenceImageType, typename TCoordRepType >
size_t
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContourType* prior ) {
	this->m_Priors.push_back( prior );
	this->m_NumberOfContours++;
    this->m_ROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentMaps.resize( this->m_NumberOfContours+1 );

	return this->m_NumberOfContours-1;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeDerivative() {
	size_t cpid = 0;
	NormalFilterPointer normalsFilter;
	SampleType sample;


	VectorType zerov; zerov.Fill(0.0);
	this->m_Derivative->FillBuffer( zerov );

	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++) {
		GradientSample s;
		PointValueType gradSum;
		PointValueType scaler = 1.0;
		sample.clear();
		// Compute mesh of normals
		normalsFilter = NormalFilterType::New();
		normalsFilter->SetInput( this->m_CurrentContours[contid] );
		normalsFilter->Update();
		ContourPointer normals = normalsFilter->GetOutput();

		typename ContourType::PointsContainerConstIterator c_it = normals->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

		PointValueType gradient;
		PointType  ci_prime;
		VectorType ni;
		typename ContourType::PointIdentifier pid;
		size_t outer_contid;

		// for every node in the mesh
		while (c_it!=c_end) {

			pid = c_it.Index();
			outer_contid = this->m_OuterList[contid][pid];

			if ( contid != outer_contid ) {
				ci_prime = c_it.Value();
				gradient =  this->GetEnergyAtPoint( ci_prime, outer_contid ) - this->GetEnergyAtPoint( ci_prime, contid );
				assert( !std::isnan(gradient) );
			} else {
				ni = zerov;
				gradient = 0.0;
			}
			s.grad = gradient;
			s.id = pid;
			sample.push_back( s );
			++c_it;
		}

		std::sort(sample.begin(), sample.end(), by_grad() );

		size_t q1 = floor( (sample.size()-1)*0.15 );
		size_t q2 = round( (sample.size()-1)*0.50 );
		size_t q3 = ceil ( (sample.size()-1)*0.85 );

		PointValueType median= sample[q2].grad;
		PointValueType quart1= sample[q1].grad;
		PointValueType quart2= sample[q3].grad;

		PointValueType maxq = ( fabs(quart1)>fabs(quart2) )?fabs(quart1):fabs(quart2);

		if( maxq >= 5.0 ) {
			scaler = 5.0 / maxq;
		}

		vnl_random rnd = vnl_random();
		for( size_t i = 0; i<q1; i++ ){
			ni = zerov;
			gradient = sample[rnd.lrand32(q1,q2)].grad * scaler;
			pid = sample[i].id;

			if ( gradient > 0.0 ) {
				// project to normal, updating transform
				normals->GetPointData( pid, &ni );         // Normal ni in point c'_i
				ni*= gradient;
			}
			normals->SetPointData( pid, ni );
			this->m_FieldInterpolator->SetPointData(cpid, ni);
			cpid++;

			sample[i].grad = gradient;
			gradSum+= gradient;
		}

		for( size_t i = q1; i<q3; i++ ){
			ni = zerov;
			gradient = sample[i].grad * scaler;
			pid = sample[i].id;

			if ( gradient > 1.0e-5 ) {
				// project to normal, updating transform
				normals->GetPointData( pid, &ni );         // Normal ni in point c'_i
				ni*= gradient;
			}
			normals->SetPointData( pid, ni );
			this->m_FieldInterpolator->SetPointData(cpid, ni);
			cpid++;

			sample[i].grad = gradient;
			gradSum+= gradient;
		}

		for( size_t i = q3; i<sample.size(); i++ ){
			ni = zerov;
			gradient = sample[rnd.lrand32(q2,q3)].grad * scaler;
			pid = sample[i].id;

			if ( gradient > 0.0 ) {
				// project to normal, updating transform
				normals->GetPointData( pid, &ni );         // Normal ni in point c'_i
				ni*= gradient;
			}
			normals->SetPointData( pid, ni );
			this->m_FieldInterpolator->SetPointData(cpid, ni);
			cpid++;

			sample[i].grad = gradient;
			gradSum+= gradient;
		}

#ifndef NDEBUG
		PointValueType minGradient = (*( sample.begin() )).grad;
		PointValueType maxGradient = (*( sample.end()-1 )).grad;
		PointValueType average = gradSum / sample.size();
		std::sort(sample.begin(), sample.end(), by_grad() );
		median= sample[q2].grad;
		quart1= sample[q1].grad;
		quart2= sample[q3].grad;
		minGradient = (*( sample.begin() )).grad;
		maxGradient = (*( sample.end()-1 )).grad;
		average = gradSum / sample.size();
		std::cout << "Grad["<< contid << "]: avg=" << average << ", max=" << maxGradient << ", min=" << minGradient << ", q1=" << quart1 << ", q2=" << quart2 << ", med=" << median << "." << std::endl;
#endif
	}

	// Interpolate sparse velocity field to targetDeformation
	this->m_FieldInterpolator->ComputeWeights();

	VectorType* gmBuffer = this->m_Derivative->GetBufferPointer();
	VectorType v;

	for( size_t gpid = 0; gpid<this->m_NumberOfNodes; gpid++ ) {
		v = this->m_FieldInterpolator->GetNodeWeight(gpid);
		if ( v.GetNorm() > 1.0e-4 ) {
			*( gmBuffer + gpid ) = v; // * (1.0/this->m_NumberOfPoints);
		}
	}

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
		if ( v.GetNorm()> 1.0e-5 ) {
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
		typename ContourType::PointsContainerConstIterator p_it = this->m_Priors[contid]->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator p_end = this->m_Priors[contid]->GetPoints()->End();

		ContourPointType ci, ci_prime;
		VectorType desp;
		size_t pid;

#ifndef NDEBUG
	double maxDesp = 0.0;
	double sumDesp = 0.0;
	size_t position =-1;
#endif

		// For all the points in the mesh
		while ( p_it != p_end ) {
			ci = p_it.Value();
			pid = p_it.Index();
			// Interpolate the value of the field in the point
			desp = this->m_FieldInterpolator->GetPointData( gpid );
			norm = desp.GetNorm();
			// Add vector to the point
			if( norm > 1.0e-3 ) {
				ci_prime = ci + desp;

				if( ! this->IsInside(ci_prime,point_idx) ) {
					itkExceptionMacro( << "Contour is outside image regions after update.\n" <<
						"\tMoving vertex [" << contid << "," << pid << "] from " << ci << " to " << ci_prime << " norm="  << desp.GetNorm() << "mm.\n");
				}

				this->m_CurrentContours[contid]->SetPoint( pid, ci_prime );

#ifndef NDEBUG
				sumDesp+=desp.GetNorm();
				if ( desp.GetNorm() > maxDesp ) {
					maxDesp = desp.GetNorm();
					position = pid;
				}
#endif
				changed++;
				gpid++;
			}
			++p_it;
		}

#ifndef NDEBUG
		std::cout << "Disp[" << contid << "]: avg="<< ( sumDesp/ (pid+1) ) << "; max=" << maxDesp << "." << std::endl;
#endif
	}


	this->m_RegionsUpdated = (changed==0);
	this->m_EnergyUpdated = (changed==0);
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	if ( !this->m_EnergyUpdated ) {
		this->m_Value = 0.0;

		double normalizer = 1.0;

		for(size_t i = 0; i<Dimension; i++)
			normalizer *= this->GetCurrentMap(0)->GetSpacing()[i];

		for( size_t roi = 0; roi < m_ROIs.size(); roi++ ) {
			ProbabilityMapConstPointer roipm = this->GetCurrentMap( roi );

	#ifndef NDEBUG
			ProbabilityMapPointer tmpmap = ProbabilityMapType::New();
			tmpmap->SetOrigin( roipm->GetOrigin() );
			tmpmap->SetRegions( roipm->GetLargestPossibleRegion() );
			tmpmap->SetDirection( roipm->GetDirection() );
			tmpmap->SetSpacing( roipm->GetSpacing() );
			tmpmap->Allocate();
			tmpmap->FillBuffer( 0.0 );

			typename ProbabilityMapType::PixelType* tmpBuffer = tmpmap->GetBufferPointer();
	#endif

			const typename ProbabilityMapType::PixelType* roiBuffer = roipm->GetBufferPointer();
			const ReferencePixelType* refBuffer = this->m_ReferenceImage->GetBufferPointer();

			size_t nPix = roipm->GetLargestPossibleRegion().GetNumberOfPixels();
			ReferencePointType pos;
			ReferencePixelType val;
			typename ProbabilityMapType::PixelType w;

			for( size_t i = 0; i < nPix; i++) {
				w = *( roiBuffer + i );
				if ( w > 0.0 ) {
					val = *(refBuffer+i);
					this->m_Value +=  w * this->GetEnergyOfSample( val, roi );
	#ifndef NDEBUG
					*(tmpBuffer+i) = val[0];
	#endif
				}
			}

	#ifndef NDEBUG
			typedef typename itk::ImageFileWriter< ProbabilityMapType > W;
			typename W::Pointer writer = W::New();
			writer->SetInput(tmpmap);
			std::stringstream ss;
			ss << "region_energy_" << roi << ".nii.gz";
			writer->SetFileName( ss.str().c_str() );
			writer->Update();
	#endif
		}
		this->m_Value = normalizer*this->m_Value;
		this->m_EnergyUpdated = true;
	}
	return this->m_Value;
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
typename FunctionalBase<TReferenceImageType, TCoordRepType>::ROIConstPointer
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetCurrentRegion( size_t idx ) {
	if(!this->m_RegionsUpdated )
		this->ComputeCurrentRegions();

	return this->m_CurrentROIs[idx];
}

template< typename TReferenceImageType, typename TCoordRepType >
const typename FunctionalBase<TReferenceImageType, TCoordRepType>::ProbabilityMapType*
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetCurrentMap( size_t idx ) {
	if(!this->m_RegionsUpdated ) {
		this->ComputeCurrentRegions();
	}

	if( this->m_CurrentMaps[idx].IsNull() ) {
		this->m_CurrentMaps[idx] = ProbabilityMapType::New();
		this->m_CurrentMaps[idx]->SetRegions(   this->m_ReferenceImage->GetLargestPossibleRegion().GetSize() );
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
::InitializeSamplingGrid() {
	typedef itk::ContinuousIndex< double, Dimension > ContinuousIndex;

	typename FieldType::SizeType size = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
	typename FieldType::SizeType exp_size;

	ContinuousIndex f_idx, l_idx;
	f_idx.Fill( -0.5 );

	for (size_t i = 0; i<Dimension; i++ ){
		l_idx[i] = size[i]-0.5;
		exp_size[i] = (unsigned int) (size[i]*this->m_SamplingFactor);
	}

	PointType first, last, origin, step;
	typename FieldType::SpacingType spacing;
	this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( f_idx, first );
	this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( l_idx, last );

	for (size_t i = 0; i<Dimension; i++ ){
		step[i] = (last[i]-first[i])/(1.0*exp_size[i]);
		spacing[i]= fabs( step[i] );
		origin[i] = first[i] + 0.5 * step[i];
	}

	this->m_ReferenceSamplingGrid = FieldType::New();
	this->m_ReferenceSamplingGrid->SetOrigin( origin );
	this->m_ReferenceSamplingGrid->SetDirection( this->m_ReferenceImage->GetDirection() );
	this->m_ReferenceSamplingGrid->SetRegions( exp_size );
	this->m_ReferenceSamplingGrid->SetSpacing( spacing );
	this->m_ReferenceSamplingGrid->Allocate();

	this->m_CurrentRegions = ROIType::New();
	this->m_CurrentRegions->SetSpacing(   this->m_ReferenceSamplingGrid->GetSpacing() );
	this->m_CurrentRegions->SetDirection( this->m_ReferenceSamplingGrid->GetDirection() );
	this->m_CurrentRegions->SetOrigin(    this->m_ReferenceSamplingGrid->GetOrigin() );
	this->m_CurrentRegions->SetRegions(   this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
	this->m_CurrentRegions->Allocate();
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


// Reorient contours to image direction in order that allowing pixel-wise computations
// ReorientFilter computes the new extent of the image if the directions
// matrix is identity. This is necessary to be able to binarize the contours
// (that are given in physical coordinates).
// See https://github.com/oesteban/ACWE-Registration/issues/92
template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeCurrentContours() {
	// Set number of control points in the sparse-dense interpolator
	this->m_NumberOfPoints = 0;
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		this->m_NumberOfPoints+= this->m_Priors[contid]->GetNumberOfPoints();
	}
	this->m_FieldInterpolator->SetN(this->m_NumberOfPoints);

	// Copy contours
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		ContourCopyPointer copy = ContourCopyType::New();
		copy->SetInput( this->m_Priors[contid] );
		copy->Update();
		this->m_CurrentContours.push_back( copy->GetOutput() );
	}

	// Cache image properties
	this->m_Origin    = this->m_ReferenceImage->GetOrigin();
	this->m_Direction = this->m_ReferenceImage->GetDirection();

	//typename ReferenceImageType::IndexType endIdx;
	//for ( size_t d = 0; d<Dimension; d++)
	//	endIdx[d] = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize()[d] -1;
	//this->m_ReferenceImage->TransformIndexToPhysicalPoint( endIdx, this->m_End );
	//DirectionType idDir; idDir.SetIdentity();
    //
	//if ( this->m_Direction != idDir ) {
	//	typedef itk::OrientImageFilter< ReferenceImageType, ReferenceImageType >   ReorientFilterType;
	//	typename ReorientFilterType::Pointer reorient = ReorientFilterType::New();
	//	reorient->UseImageDirectionOn();
	//	reorient->SetDesiredCoordinateDirection(idDir);
	//	reorient->SetInput( this->m_ReferenceImage );
	//	reorient->Update();
	//	ReferenceImagePointer reoriented = reorient->GetOutput();
	//	ReferencePointType newOrigin = reoriented->GetOrigin();
    //
	//	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
	//		ContourPointer c = this->m_CurrentContours[contid];
	//		PointsIterator c_it  = c->GetPoints()->Begin();
	//		PointsIterator c_end = c->GetPoints()->End();
	//		ContourPointType ci, ci_new;
	//		size_t pid;
	//		while( c_it != c_end ) {
	//			pid = c_it.Index();
	//			ci = c_it.Value();
	//			ci_new = (this->m_Direction* ( ci - this->m_Origin.GetVectorFromOrigin() )) + newOrigin.GetVectorFromOrigin();
	//			c->SetPoint( pid, ci_new );
	//			++c_it;
	//		}
	//	}
	//}

	// Fill in interpolator points
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		PointsIterator c_it  = this->m_CurrentContours[contid]->GetPoints()->Begin();
		PointsIterator c_end = this->m_CurrentContours[contid]->GetPoints()->End();
		PointType ci;
		size_t pid;
		size_t cpid = 0;
		while( c_it != c_end ) {
			pid = c_it.Index();
			ci = c_it.Value();
			this->m_FieldInterpolator->SetPoint( cpid, ci );
			++c_it;
			cpid++;
		}
	}
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
			meshFilter->SetInput(     this->m_CurrentContours[idx]);
			meshFilter->Update();
			tempROI = meshFilter->GetOutput();
		} else {
			tempROI = ROIType::New();
			tempROI->SetSpacing(   this->m_ReferenceSamplingGrid->GetSpacing() );
			tempROI->SetDirection( this->m_ReferenceSamplingGrid->GetDirection() );
			tempROI->SetOrigin(    this->m_ReferenceSamplingGrid->GetOrigin() );
			tempROI->SetRegions(   this->m_ReferenceSamplingGrid->GetLargestPossibleRegion().GetSize() );
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

	this->m_RegionsUpdated = true;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeOuterRegions() {
	ContourOuterRegions outerVect;

	if( this->m_NumberOfContours>1 ) {
		// Set up ROI interpolator
		typename ROIInterpolatorType::Pointer interp = ROIInterpolatorType::New();
		interp->SetInputImage( this->m_CurrentRegions );

		// Set up outer regions
		for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
			// Compute mesh of normals
			NormalFilterPointer normalsFilter = NormalFilterType::New();
			normalsFilter->SetInput( this->m_Priors[contid] );
			normalsFilter->Update();
			ContourPointer normals = normalsFilter->GetOutput();
			outerVect.resize( normals->GetNumberOfPoints() );

			typename ContourType::PointsContainerConstIterator c_it  = normals->GetPoints()->Begin();
			typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

			ContourPointType ci;
			VectorType v;
			VectorType ni;

			size_t pid;
			while( c_it != c_end ) {
				pid = c_it.Index();
				ci = c_it.Value();
				//ci = normals->GetPoint(pid);
				normals->GetPointData( pid, &ni );
				outerVect[pid] =  interp->Evaluate( ci - ni * 0.5 );
				++c_it;
			}
			this->m_OuterList.push_back( outerVect );
		}
	} else {
		outerVect.resize(this->m_CurrentContours[0]->GetNumberOfPoints());
		std::fill( outerVect.begin(), outerVect.end(), 1 );
		this->m_OuterList.push_back( outerVect );
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeInterpolatorGrid() {
	if( this->m_Derivative.IsNotNull() ) {
		PointType uk;
		this->m_NumberOfNodes = this->m_Derivative->GetLargestPossibleRegion().GetNumberOfPixels();
		this->m_FieldInterpolator->SetNumberOfParameters( this->m_NumberOfNodes );

		for ( size_t gid = 0; gid < this->m_NumberOfNodes; gid++ ) {
			this->m_Derivative->TransformIndexToPhysicalPoint( this->m_Derivative->ComputeIndex( gid ), uk );
			this->m_FieldInterpolator->SetNode( gid, uk );
		}
	} else {
		itkWarningMacro( << "No parametrization (deformation field grid) was defined.");
	}

	typename FieldType::SpacingType sigma = this->m_Derivative->GetSpacing();
	this->m_FieldInterpolator->SetSigma( sigma );
}

}

#endif /* LEVELSETSBASE_HXX_ */
