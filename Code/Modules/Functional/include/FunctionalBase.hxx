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

#ifndef FUNCTIONALBASE_HXX_
#define FUNCTIONALBASE_HXX_

#include "FunctionalBase.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <numeric>
#include <vnl/vnl_random.h>
#include <itkImageAlgorithm.h>
#include <itkOrientImageFilter.h>
#include <itkContinuousIndex.h>

#include <itkMeshFileWriter.h>

#define MAX_GRADIENT 20.0
#define MIN_GRADIENT 1.0e-8

namespace rstk {

template< typename TReferenceImageType, typename TCoordRepType >
FunctionalBase<TReferenceImageType, TCoordRepType>
::FunctionalBase():
 m_NumberOfContours(0),
 m_NumberOfRegions(1),
 m_NumberOfVertices(0),
 m_SamplingFactor(2),
 m_DecileThreshold(0.05),
 m_DisplacementsUpdated(true),
 m_EnergyUpdated(false),
 m_RegionsUpdated(false),
 m_ApplySmoothing(false),
 m_UseBackground(false),
 m_Value(0.0),
 m_MaxEnergy(0.0),
 m_LastROI(0)
 {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_Sigma.Fill(0.0);
	this->m_Interp = InterpolatorType::New();
	this->m_MaskInterp = MaskInterpolatorType::New();
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::Initialize() {
	this->ParseSettings();

	if( this->m_ApplySmoothing ) {
		SmoothingFilterPointer s = SmoothingFilterType::New();
		s->SetInput( this->m_ReferenceImage );
		s->SetSigmaArray( this->m_Sigma );
		s->Update();
		this->m_ReferenceImage = s->GetOutput();
	}

	if (this->m_BackgroundMask.IsNull()) {
		this->m_UseBackground = false;
		ProbabilityMapPointer p = ProbabilityMapType::New();
		p->SetRegions(this->m_ReferenceImage->GetLargestPossibleRegion());
		p->SetSpacing(this->m_ReferenceImage->GetSpacing());
		p->SetOrigin(this->m_ReferenceImage->GetOrigin());
		p->SetDirection(this->m_ReferenceImage->GetDirection());
		p->Allocate();
		p->FillBuffer(0.0);
		this->m_BackgroundMask = p;
	}

	// Initialize corresponding ROI /////////////////////////////
	// Check that high-res reference sampling grid has been initialized
	if ( this->m_ReferenceSamplingGrid.IsNull() ) {
			this->InitializeSamplingGrid();
	}

	// Initialize interpolators
	this->m_Interp->SetInputImage( this->m_ReferenceImage );
	this->m_MaskInterp->SetInputImage(this->m_BackgroundMask);

	// Compute and set regions in m_ROIs
	this->ComputeCurrentRegions();
	for( size_t id = 0; id < m_ROIs.size(); id++) {
		this->m_ROIs[id] = this->m_CurrentROIs[id];
	}
	this->m_LastROI = m_ROIs.size() - 1;

	// Compute the outer region in each vertex
	this->InitializeContours();

	this->m_RegionValue.SetSize(this->m_NumberOfRegions);
	this->m_RegionValue.Fill(itk::NumericTraits<MeasureType>::infinity());

	//this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
}

template< typename TReferenceImageType, typename TCoordRepType >
size_t
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::VectorContourType* prior ) {
	this->m_Offsets.push_back( this->m_NumberOfVertices );
	this->m_Priors.push_back( prior );

	ContourCopyPointer copy = ContourCopyType::New();
	copy->SetInput( prior );
	copy->Update();
	this->m_CurrentContours.push_back( copy->GetOutput() );

	// Increase number of off-grid nodes to set into the sparse-dense interpolator
	this->m_NumberOfVertices+= prior->GetNumberOfPoints();

	this->m_NormalsFilter.push_back(NormalFilterType::New());
	this->m_NormalsFilter[this->m_NumberOfContours]->SetWeight(NormalFilterType::AREA);
	this->m_NormalsFilter[this->m_NumberOfContours]->SetInput( this->m_CurrentContours[this->m_NumberOfContours] );
	this->m_NormalsFilter[this->m_NumberOfContours]->Update();

	ContourCopyPointer copygrad = ContourCopyType::New();
	copygrad->SetInput( this->m_NormalsFilter[this->m_NumberOfContours]->GetOutput() );
	copygrad->Update();
	this->m_Gradients.push_back( copygrad->GetOutput() );

	this->m_NumberOfContours++;
	this->m_NumberOfRegions++;
    this->m_ROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentMaps.resize( this->m_NumberOfContours+1 );
	return this->m_NumberOfContours-1;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeDerivative(PointValueType* grad, ScalesType scales) {
	this->UpdateContour();
	size_t nvertices = this->m_ValidVertices.size();
	size_t fullsize = nvertices * Dimension;
	PointIdentifier pid;      // universal id of vertex
	PointIdentifier cpid;     // id of vertex in its contour

	ContourPointType ci_prime;
	PointValueType gi[nvertices];
	PointValueType sum = 0.0;
	std::vector< PointValueType > sample;
	ROIPixelType ocid;
	ROIPixelType icid = 0;
	double wi = 0.0;

	std::vector< NormalFilterAreasContainer > areas;
	std::vector< ContourPointer > normals;
	for (size_t i = 0; i < this->m_NumberOfContours; i++ ) {
		areas.push_back(this->m_NormalsFilter[i]->GetVertexAreaContainer());
		normals.push_back (this->m_NormalsFilter[i]->GetOutput() );
	}

	for(size_t vvid = 0; vvid < nvertices; vvid++ ) {
		icid = this->m_InnerRegion[vvid];
		ocid = this->m_OuterRegion[vvid];
		pid = this->m_ValidVertices[vvid];
		cpid = pid - this->m_Offsets[icid];
		this->m_CurrentContours[icid]->GetPoint(cpid, &ci_prime); // Get c'_i
		wi = areas[icid][cpid];
		gi[vvid] = this->EvaluateGradient( ci_prime, ocid, icid )  * wi;
		sample.push_back(gi[vvid]);
		sum+= fabs(gi[vvid]);
	}

	std::sort(sample.begin(), sample.end());

	this->m_GradientStatistics[0] = sample.front();
	this->m_GradientStatistics[1] = sample[int(0.05 * (sample.size()-1))];
	this->m_GradientStatistics[2] = sample[int(0.25 * (sample.size()-1))];
	this->m_GradientStatistics[3] = sample[int(0.50 * (sample.size()-1))];
	this->m_GradientStatistics[4] = sample[int(0.75 * (sample.size()-1))];
	this->m_GradientStatistics[5] = sample[int(0.95 * (sample.size()-1))];
	this->m_GradientStatistics[6] = sample.back();

	VectorType ni, v;
	icid = 0;
	PointValueType g;
	for(size_t vvid = 0; vvid < nvertices; vvid++ ) {
		pid = this->m_ValidVertices[vvid];
		if( pid == this->m_Offsets[icid + 1] ) {
			icid++;
		}

		cpid = pid - this->m_Offsets[icid];
		ni.Fill(0.0);
		normals[icid]->GetPointData(cpid, &ni);

		g = gi[vvid];
		if ( g > this->m_GradientStatistics[5] ) g = this->m_GradientStatistics[5];
		if ( g < this->m_GradientStatistics[1] ) g = this->m_GradientStatistics[1];

		for( size_t i = 0; i < Dimension; i++ ) {
			v[i] = scales[i] * g * ni[i];
			grad[vvid + i * nvertices] = static_cast<float>(v[i]);
		}

		this->m_Gradients[icid]->GetPointData()->SetElement( cpid, v );
	}

}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::UpdateContour() {
	if( this->m_DisplacementsUpdated ) {
		return;
	}

	MeasureType norm;
	ContinuousIndex point_idx;
	size_t changed = 0;
	size_t gpid = 0;
	std::vector< size_t > invalid;
	VectorType zerov; zerov.Fill(0.0);


	std::fill(this->m_OffMaskVertices.begin(), this->m_OffMaskVertices.end(), 0);

	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++ ) {
		typename VectorContourType::PointsContainerConstIterator p_it = this->m_Priors[contid]->GetPoints()->Begin();
		typename VectorContourType::PointsContainerConstIterator p_end = this->m_Priors[contid]->GetPoints()->End();
		PointsContainerPointer curPoints = this->m_CurrentContours[contid]->GetPoints();
		PointsContainerPointer curGradPoints = this->m_Gradients[contid]->GetPoints();
		PointDataContainerPointer gradData = this->m_Gradients[contid]->GetPointData();

		ContourPointType ci, ci_prime;
		VectorType disp;
		size_t pid;

		// For all the points in the mesh
		while ( p_it != p_end ) {
			ci = p_it.Value();
			pid = p_it.Index();
			disp = this->m_CurrentDisplacements->GetElement( gpid ); // Get the interpolated value of the field in the point
			norm = disp.GetNorm();

			if( norm > 1.0e-8 ) {
				ci_prime = ci + disp; // Add displacement vector to the point
				if( ! this->CheckExtent(ci_prime,point_idx) ) {
					invalid.push_back( gpid );
					this->InvokeEvent( WarningEvent() );
				}
				curPoints->SetElement( pid, ci_prime );
				curGradPoints->SetElement( pid, ci_prime );
				changed++;
			}

			if ( (1.0 - this->m_MaskInterp->Evaluate(ci_prime)) < 1.0e-5 ) {
				this->m_OffMaskVertices[contid]++;
			}

			this->m_Gradients[contid]->GetPointData()->SetElement( pid, zerov );
			++p_it;
			gpid++;
		}
	}

	if ( invalid.size() > 0 ) {
		itkWarningMacro(<< "a total of " << invalid.size() << " mesh nodes were to be moved off the image domain." );
	}

	UpdateNormals();

	this->m_DisplacementsUpdated = true;
	this->m_RegionsUpdated = (changed==0);
	this->m_EnergyUpdated = (changed==0);
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::UpdateNormals() {
	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++ ) {
		// Compute mesh of normals
		this->m_NormalsFilter[contid]->Update();
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	if ( !this->m_EnergyUpdated ) {
		this->m_Value = 0.0;
		this->m_RegionValue.Fill(0.0);

		double vxvol = 1.0;
		for(size_t i = 0; i<Dimension; i++)
			vxvol *= this->GetCurrentMap(0)->GetSpacing()[i];

		size_t nPix = this->GetCurrentMap(0)->GetLargestPossibleRegion().GetNumberOfPixels();
		size_t nrois = m_ROIs.size();
		size_t lastroi = nrois;

		const ReferencePixelType* refBuffer = this->m_ReferenceImage->GetBufferPointer();
		const typename ProbabilityMapType::PixelType* bgBuffer = this->m_BackgroundMask->GetBufferPointer();
		const typename ProbabilityMapType::PixelType* roiBuffer[nrois];

		for( size_t roi = 0; roi < lastroi; roi++ ) {
			roiBuffer[roi] = this->GetCurrentMap(roi)->GetBufferPointer();
		}

		ReferencePointType pos;
		ReferencePixelType val;
		typename ProbabilityMapType::PixelType w;
		typename ProbabilityMapType::PixelType bgw;
		double regionVol[nrois];
		regionVol[lastroi] = 0.0;

		MeasureType e;
		for( size_t i = 0; i < nPix; i++) {
			bgw = *(bgBuffer + i);

			if (bgw < 1.0) {
				val = *(refBuffer+i);
				for( size_t roi = 0; roi < lastroi; roi++ ) {
					w = *( roiBuffer[roi] + i );
					if ( w < 1.0e-8 ) {
						continue;
					}

					if (bgw > 1.0e-3 && roi < (nrois - 1)) {
						e = this->m_MaxEnergy;
					} else {
						e = this->GetEnergyOfSample( val, roi );
					}
					this->m_RegionValue[roi]+= w * vxvol * e;
					regionVol[roi]+= w * vxvol;
				}
			}
		}

		this->m_Value = 0.0;
		for( size_t roi = 0; roi < nrois; roi++ ) {
			this->m_RegionValue[roi]+= regionVol[roi] * this->GetEnergyOffset(roi);
			this->m_Value+= this->m_RegionValue[roi];
		}

		this->m_EnergyUpdated = true;
	}
	return this->m_Value;
}

template< typename TReferenceImageType, typename TCoordRepType >
inline bool
FunctionalBase<TReferenceImageType, TCoordRepType>
::CheckExtent( typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContourPointType& p, typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContinuousIndex& idx) const {
	ReferencePointType ref;
	ref.CastFrom ( p );
	bool isInside = this->m_ReferenceImage->TransformPhysicalPointToContinuousIndex( ref , idx );

	if(!isInside) {
		for ( size_t i = 0; i<Dimension; i++) {
			if ( idx[i] < 0.0 ) {
				p.SetElement(i, this->m_FirstPixelCenter[i] );
			}
			else if ( idx[i] > (this->m_ReferenceSize[i] -1) ) {
				p.SetElement(i, this->m_LastPixelCenter[i] );
			}
		}
	}

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
		this->m_CurrentMaps[idx]->SetRegions(   this->m_ReferenceSize );
		this->m_CurrentMaps[idx]->SetOrigin(    this->m_FirstPixelCenter );
		this->m_CurrentMaps[idx]->SetDirection( this->m_Direction );
		this->m_CurrentMaps[idx]->SetSpacing(   this->m_ReferenceSpacing );
		this->m_CurrentMaps[idx]->Allocate();
		this->m_CurrentMaps[idx]->FillBuffer( 0.0 );
	}

	// Resample to reference image resolution
	ResampleROIFilterPointer resampleFilter = ResampleROIFilterType::New();
	resampleFilter->SetInput( this->m_CurrentROIs[idx] );
	resampleFilter->SetSize( this->m_ReferenceSize );
	resampleFilter->SetOutputOrigin(    this->m_FirstPixelCenter );
	resampleFilter->SetOutputSpacing(   this->m_ReferenceSpacing );
	resampleFilter->SetOutputDirection( this->m_Direction );
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
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeSamplingGrid() {
	typename FieldType::SizeType exp_size;

	for (size_t i = 0; i<Dimension; i++ ){
		exp_size[i] = (unsigned int) (this->m_ReferenceSize[i]*this->m_SamplingFactor);
	}

	PointType firstPixelCenter;
	VectorType step;
	typename FieldType::SpacingType spacing;

	for (size_t i = 0; i<Dimension; i++ ){
		step[i] = (this->m_End[i] - this->m_Origin[i]) / (1.0*exp_size[i]);
		spacing[i]= fabs( step[i] );
		firstPixelCenter[i] = this->m_Origin[i] + 0.5 * step[i];
	}

	this->m_ReferenceSamplingGrid = FieldType::New();
	this->m_ReferenceSamplingGrid->SetOrigin( firstPixelCenter );
	this->m_ReferenceSamplingGrid->SetDirection( this->m_Direction );
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

		double total = 0.0;
		for( size_t pix = 0; pix < nPix; pix++ ) {
			if( *(regionsBuffer+pix) == unassigned && *( roiBuffer + pix )==1 ) {
				*(regionsBuffer+pix) = idx;
				total += 1;
			} else {
				*( roiBuffer + pix ) = 0;
			}
		}

		if (total == 0.0) {
			itkWarningMacro(<< " ROI " << idx << " is empty.")
		}

		this->m_CurrentROIs[idx] = tempROI;
	}
	this->m_RegionsUpdated = true;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeContours() {
	this->m_OffMaskVertices.resize(this->m_NumberOfContours);
	std::fill(this->m_OffMaskVertices.begin(), this->m_OffMaskVertices.end(), 0);

	VectorType zerov; zerov.Fill(0.0);
	this->m_CurrentDisplacements = PointDataContainer::New();
	this->m_CurrentDisplacements->Reserve( this->m_NumberOfVertices );

	if( this->m_NumberOfRegions > 2 ) {
		// Set up ROI interpolator
		typename ROIInterpolatorType::Pointer interp = ROIInterpolatorType::New();
		interp->SetInputImage( this->m_CurrentRegions );

		PointIdentifier tpid = 0;

		// Set up outer regions
		for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
			// Compute mesh of normals
			NormalFilterPointer normalsFilter = NormalFilterType::New();
			normalsFilter->SetInput( this->m_CurrentContours[contid] );
			normalsFilter->Update();
			ContourPointer normals = normalsFilter->GetOutput();

			PointsConstIterator c_it  = normals->GetPoints()->Begin();
			PointsConstIterator c_end = normals->GetPoints()->End();

			PointType ci;
			VectorType v;
			VectorType ni;
			size_t pid;
			float step = 1;
			while( c_it != c_end ) {
				pid = c_it.Index();
				ci = c_it.Value();

				this->m_Vertices.push_back( ci );
				this->m_CurrentDisplacements->SetElement(tpid, zerov);

				if ( (1.0 - this->m_MaskInterp->Evaluate(ci)) < 1.0e-5 ) {
					this->m_OffMaskVertices[contid]++;
				}


				this->m_Gradients[contid]->GetPointData()->SetElement( pid, zerov );
				normals->GetPointData( pid, &ni );
				ROIPixelType inner = interp->Evaluate( ci + ni );
				ROIPixelType outer = interp->Evaluate( ci - ni );

				if ((inner == contid)) {
					while (outer==inner) {
						outer = interp->Evaluate( ci - (ni * step * 0.1) );
						if (step == 9)
							break;
						step++;
					}

					if(outer!=inner) {
						this->m_ValidVertices.push_back(tpid);
						this->m_OuterRegion.push_back(outer);
						this->m_InnerRegion.push_back(inner);
					}
				}

				++c_it;
				tpid++;
			}
		}
	} else {
		this->m_InnerRegion.resize(this->m_NumberOfVertices);
		std::fill(this->m_InnerRegion.begin(), this->m_InnerRegion.end(), 0);
		this->m_OuterRegion.resize(this->m_NumberOfVertices);
		std::fill(this->m_OuterRegion.begin(), this->m_OuterRegion.end(), 1);
	}
	std::cout << "Valid vertices: " << this->m_ValidVertices.size() << " of " << this->m_Vertices.size() << "." << std::endl;
}

template< typename TReferenceImageType, typename TCoordRepType >
double
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputePointArea(const PointIdentifier & iId, VectorContourType *mesh ) {
	QEType* edge = mesh->FindEdge( iId );
	QEType* temp = edge;
	CellIdentifier cell_id(0);
	double totalArea = 0.0;
	ContourPointType pt[3];
	typedef typename PolygonType::PointIdIterator PolygonPointIterator;

	do {
		cell_id = temp->GetLeft();

		if ( cell_id != VectorContourType::m_NoFace ) {
			PolygonType *poly = dynamic_cast< PolygonType * >(
					mesh->GetCells()->GetElement(cell_id) );
			PolygonPointIterator pit = poly->PointIdsBegin();

			for(size_t k = 0; pit!= poly->PointIdsEnd(); ++pit, k++ ) {
				pt[k] = mesh->GetPoint( *pit );
			}

			totalArea += TriangleType::ComputeArea(pt[0], pt[1], pt[2]);
		}

		temp = temp->GetOnext();
	} while ( temp != edge );
	return fabs(totalArea * 0.33);
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddOptions( SettingsDesc& opts ) {
	opts.add_options()
			("smoothing", bpo::value< float > (), "apply isotropic smoothing filter on target image, with kernel sigma=S mm.")
			("smooth-auto", bpo::bool_switch(), "apply isotropic smoothing filter on target image, with automatic computation of kernel sigma.")
			//("use-background", bpo::bool_switch(), "consider last ROI as background and do not compute descriptors.")
			("decile-threshold,d", bpo::value< float > (), "set (decile) threshold to consider a computed gradient as outlier (ranges 0.0-0.5)");
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ParseSettings() {
	if( this->m_Settings.count( "smoothing" ) ) {
		this->m_ApplySmoothing = true;
		bpo::variable_value v = this->m_Settings["smoothing"];
		this->m_Sigma.Fill( v.as<float> () );
	}

	if( this->m_Settings.count( "smooth-auto" ) ) {
		bpo::variable_value v = this->m_Settings["smooth-auto"];
		if ( v.as<bool>() ) {
			this->m_ApplySmoothing= v.as<bool> ();
			this->m_Sigma.Fill( 0.0 );
		}
	}

	//if( this->m_Settings.count( "use-background" ) ) {
	//	bpo::variable_value v = this->m_Settings["use-background"];
	//	if ( v.as<bool>() ) {
	//		this->SetUseBackground(true);
	//	}
	//}

	if( this->m_Settings.count( "decile-threshold") ) {
		bpo::variable_value v = this->m_Settings["decile-threshold"];
		this->SetDecileThreshold( v.as<float> () );
	}
	this->Modified();
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::SetCurrentDisplacements( const VNLVectorContainer& vals ) {
	if ( vals.Size() != Dimension ) {
		itkExceptionMacro(<< "vals contains a wrong number of columns");
	}
	size_t npoints = vals[0].size();

	if ( npoints != this->m_NumberOfVertices ) {
		itkExceptionMacro( << "vals contains a wrong number of vectors");
	}

	VectorType new_ci, old_ci;
	MeasureType norm;
	size_t modified = 0;
	for( size_t id = 0; id<npoints; id++ ) {
		new_ci.Fill( 0.0 );
		norm = 0.0;

		for( size_t d = 0; d < Dimension; d++) {
			new_ci[d] = vals[d][id];
		}

		old_ci = this->m_CurrentDisplacements->GetElement( id );
		norm = (new_ci-old_ci).GetNorm();

		if ( norm > 1.0e-8 ) {
			modified++;
			this->m_CurrentDisplacements->SetElement( id, new_ci );
		}
	}

	this->m_DisplacementsUpdated = (modified==0);

}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::SetReferenceImage ( const ReferenceImageType * _arg ) {
	itkDebugMacro("setting ReferenceImage to " << _arg);

	if ( this->m_ReferenceImage != _arg ) {
		typedef itk::OrientImageFilter< ReferenceImageType, ReferenceImageType >  Orienter;
		typename Orienter::Pointer orient = Orienter::New();
		orient->UseImageDirectionOn();
		orient->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
		orient->SetInput(_arg);
		orient->Update();
		ReferenceImagePointer ref = orient->GetOutput();
		ReferenceSizeType size = ref->GetLargestPossibleRegion().GetSize();

		// ReferenceImagePointer newref = ReferenceImageType::New();
		// newref->CopyInformation(ref);
		// newref->SetRegions(size);
		// newref->Allocate();

		// ReferenceIndexType i_old;
		// ReferenceIndexType i_new;
		// //size_t off_old, off_new;
		// for(size_t i = 0; i< size[0]; i++) {
		// 	i_old[0] = i;
		// 	i_new[0] = size[0] - i - 1;
		// 	for(size_t j = 0; j < size[1]; j++) {
		// 		i_old[1] = j;
		// 		i_new[1] = size[1] - j - 1;
		// 		for(size_t k = 0; k < size[2]; k++) {
		// 			i_old[2] = k;
		// 			i_new[2] = k;
		// 			newref->SetPixel(i_new, ref->GetPixel(i_old));
		// 		}
		// 	}
		// }


		DirectionType idmat; idmat.SetIdentity();
		DirectionType itk; itk.SetIdentity();
		itk(0,0) = -1.0; itk(1,1) = -1.0;

		PointType neworig = itk * ref->GetOrigin();
		ref->SetDirection(idmat);
		ref->SetOrigin(neworig);
		this->m_ReferenceImage = ref;




		// Cache image properties
		this->m_FirstPixelCenter  = this->m_ReferenceImage->GetOrigin();
		this->m_Direction = this->m_ReferenceImage->GetDirection();
		this->m_ReferenceSize = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
		this->m_ReferenceSpacing = this->m_ReferenceImage->GetSpacing();

		ContinuousIndex tmp_idx;
		tmp_idx.Fill( -0.5 );
		this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( tmp_idx, this->m_Origin );

		for ( size_t dim = 0; dim<FieldType::ImageDimension; dim++)  tmp_idx[dim]= this->m_ReferenceSize[dim]-1.0;
		this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( tmp_idx, this->m_LastPixelCenter );

		for ( size_t dim = 0; dim<FieldType::ImageDimension; dim++)  tmp_idx[dim]= this->m_ReferenceSize[dim]- 0.5;
		this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( tmp_idx, this->m_End );

		this->Modified();
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::SetBackgroundMask (const ProbabilityMapType * _arg) {
	if ( this->GetDebug() && ::itk::Object::GetGlobalWarningDisplay() ) {
		std::ostringstream itkmsg;
		itkmsg << "Debug: In " "/home/oesteban/workspace/ACWE-Reg/Code/Modules/Functional/include/FunctionalBase.h" ", line " << 333 << "\n" \
				<< this->GetNameOfClass() << " (" << this << "): " "setting " << "BackgroundMask" " to " << _arg       \
				<< "\n\n";
		::itk::OutputWindowDisplayDebugText( itkmsg.str().c_str() );
	}
	if ( this->m_BackgroundMask != _arg ) {
		DirectionType olddir = _arg->GetDirection();
		PointType oldorig = _arg->GetOrigin();
		ProbabilityMapPointer msk = ProbabilityMapType::New();
		msk->CopyInformation( _arg );
		msk->SetRegions( _arg->GetLargestPossibleRegion() );
		msk->Allocate();
		msk->FillBuffer(0.0);

		size_t nPix = msk->GetLargestPossibleRegion().GetNumberOfPixels();
		float* mbuffer = msk->GetBufferPointer();
		const float* sbuffer = _arg->GetBufferPointer();
		float w;

		for (size_t i = 0; i<nPix; i++) {
			w = *(sbuffer + i);

			if (w < 1.0) {
				w = 1.0 - w;
				if ( w > 1.0e-3 )
					if (w > 1.0)
						w = 1.0;

					*(mbuffer + i) = w;
			}
		}

		DirectionType ident; ident.SetIdentity();
		msk->SetDirection(ident);

		itk::Vector<double, Dimension> ovect = oldorig.GetVectorFromOrigin();
		msk->SetOrigin(olddir * ovect);

		if ( (msk->GetLargestPossibleRegion() != this->m_ReferenceImage->GetLargestPossibleRegion()) ||
				(msk->GetOrigin() != this->m_ReferenceImage->GetOrigin())) {
			ProbmapResamplePointer res = ProbmapResampleType::New();
			res->SetInput(msk);
			res->SetSize(this->m_ReferenceSize);
			res->SetOutputSpacing(this->m_ReferenceSpacing);
			res->SetOutputOrigin(this->m_FirstPixelCenter);
			res->SetOutputDirection(this->m_Direction);
			res->SetDefaultPixelValue( 0.0 );
			res->Update();

			this->m_BackgroundMask = res->GetOutput();
		}
		else {
			this->m_BackgroundMask = msk;
		}
		this->Modified();
	}
}


template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::EvaluateGradient( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point,
		size_t outer_roi, size_t inner_roi ) const {
	if ( outer_roi == inner_roi ) {
		return 0.0;
	}
	ReferencePixelType value = this->m_Interp->Evaluate( point );
	MeasureType gin  = this->GetEnergyOfSample( value, inner_roi );
	MeasureType gout = this->GetEnergyOfSample( value, outer_roi );

	float isOutside = this->m_MaskInterp->Evaluate( point );
	if (isOutside > 1.0e-3) {
		if(isOutside > 1.0) isOutside = 1.0;
		gout = isOutside * this->m_MaxEnergy;
	}

	MeasureType grad = gin - gout;
	grad = (fabs(grad)>MIN_GRADIENT)?grad:0.0;
	return grad;
}


template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetEnergyAtPoint( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point, size_t roi ) const {
	ReferencePixelType value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	if (isOutside > 1.0)
		isOutside = 1.0;
	else if (isOutside < 1.0e-3)
		isOutside = 0.0;
	float factor = 1.0 * (this->m_NumberOfRegions - roi) / this->m_NumberOfRegions;
	return this->GetEnergyOfSample( value, roi ) + isOutside * factor * this->m_MaxEnergy;
}

template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetEnergyAtPoint( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point, size_t roi,
		            typename FunctionalBase<TReferenceImageType, TCoordRepType>::ReferencePixelType & value) const {
	value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	if (isOutside > 1.0)
		isOutside = 1.0;
	else if (isOutside < 1.0e-3)
		isOutside = 0.0;
	float factor = 1.0 * (this->m_NumberOfRegions - roi) / this->m_NumberOfRegions;
	return this->GetEnergyOfSample( value, roi ) + isOutside * factor * this->m_MaxEnergy;
}

}

#endif /* FUNCTIONALBASE_HXX_ */
