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
#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

#include <itkMeshFileWriter.h>
#include <itkImageAlgorithm.h>
#include <itkOrientImageFilter.h>
#include <itkContinuousIndex.h>

#define MAX_GRADIENT 20.0
#define MIN_GRADIENT 1.0e-5

namespace rstk {


template< typename TReferenceImageType, typename TCoordRepType >
FunctionalBase<TReferenceImageType, TCoordRepType>
::FunctionalBase():
 m_NumberOfContours(0),
 m_NumberOfRegions(1),
 m_NumberOfPoints(0),
 m_NumberOfNodes(0),
 m_SamplingFactor(2),
 m_DecileThreshold(0.05),
 m_DisplacementsUpdated(true),
 m_EnergyUpdated(false),
 m_RegionsUpdated(false),
 m_ApplySmoothing(false),
 m_UseBackground(false),
 m_Value(0.0),
 m_MaxEnergy(0.0)
 {
	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_Sigma.Fill(0.0);
	this->m_Interp = InterpolatorType::New();
	this->m_MaskInterp = ProbmapInterpolatorType::New();
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

	// Initialize interpolators
	this->m_Interp->SetInputImage( this->m_ReferenceImage );
	this->m_MaskInterp->SetInputImage(this->m_BackgroundMask);
}

template< typename TReferenceImageType, typename TCoordRepType >
size_t
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContourType* prior ) {
	this->m_Priors.push_back( prior );

	// Increase number of off-grid nodes to set into the sparse-dense interpolator
	this->m_NumberOfPoints+= prior->GetNumberOfPoints();

	this->m_NumberOfContours++;
	this->m_NumberOfRegions++;
    this->m_ROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentROIs.resize( this->m_NumberOfContours+1 );
    this->m_CurrentMaps.resize( this->m_NumberOfContours+1 );
	return this->m_NumberOfContours-1;
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::VNLVectorContainer
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeDerivative() {
	size_t cpid = 0;
	NormalFilterPointer normalsFilter;
	SampleType sample;
	VectorType zerov; zerov.Fill(0.0);
	this->UpdateContour();

	VNLVectorContainer gradVector;

	for( size_t i = 0; i<Dimension; i++ ) {
		gradVector[i] = VNLVector( this->m_NumberOfPoints );
		gradVector[i].fill(0.0);
	}

	for( size_t in_cid = 0; in_cid < this->m_NumberOfContours; in_cid++) {
		sample.clear();
		double wi = 0.0;
		PointValueType totalArea = 0.0;
		PointValueType gradSum = 0.0;

		// Compute mesh of normals
		normalsFilter = NormalFilterType::New();
		normalsFilter->SetInput( this->m_CurrentContours[in_cid] );
		normalsFilter->Update();
		ContourPointer normals = normalsFilter->GetOutput();

		typename ContourType::PointsContainerConstIterator c_it = normals->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

		PointType  ci_prime;
		VectorType ni;
		PointValueType gi;
		typename ContourType::PointIdentifier pid;
		size_t out_cid;

		// for every node in the mesh: compute gradient, assign cid and gid.
		while (c_it!=c_end) {
			ni = zerov;
			gi = 0.0;
			wi = 0.0;

			pid = c_it.Index();
			out_cid = this->m_OuterList[in_cid][pid];

			if ( in_cid != out_cid ) {
				ci_prime = c_it.Value();
				normals->GetPointData( pid, &ni );           // Normal ni in point c'_i
				wi = this->ComputePointArea( pid, normals );  // Area of c'_i
				gi = this->EvaluateGradient( ci_prime, out_cid, in_cid );
				totalArea+=wi;
				gradSum+=gi;
			}
			sample.push_back( GradientSample( gi, wi, ni, pid, cpid, in_cid ) );
			++c_it;
			cpid++;
		}

		PointValueType gradient;
		ShapeGradientPointer gradmesh = this->m_Gradients[in_cid];
		gradSum = 0.0;
		for( size_t i = 0; i< sample.size(); i++) {
			if ( sample[i].w > 0.0 ) {
				gradient = sample[i].grad;
				sample[i].grad = gradient;
				gradSum+= gradient;
				ni = gradient * sample[i].normal;  // Project to normal

				for( size_t dim = 0; dim<Dimension; dim++ ) {
					if( fabs(ni[dim]) > MIN_GRADIENT )
						gradVector[dim][sample[i].gid] = ni[dim];
				}
			} else {
				sample[i].normal = zerov;
				sample[i].grad = 0.0;
				sample[i].w = 0.0;
				ni = zerov;
			}

			gradmesh->GetPointData()->SetElement( sample[i].cid, ni );
		}
	}

	return gradVector;
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

	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++ ) {
		typename ContourType::PointsContainerConstIterator p_it = this->m_Priors[contid]->GetPoints()->Begin();
		typename ContourType::PointsContainerConstIterator p_end = this->m_Priors[contid]->GetPoints()->End();
		PointsContainerPointer curPoints = this->m_CurrentContours[contid]->GetPoints();

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
				changed++;
			}
			++p_it;
			gpid++;
		}

		ShapeCopyPointer copyShape = ShapeCopyType::New();
		copyShape->SetInput( this->m_CurrentContours[contid] );
		copyShape->Update();
		this->m_Gradients[contid] = copyShape->GetOutput();
	}

	if ( invalid.size() > 0 ) {
		itkWarningMacro(<< "a total of " << invalid.size() << " mesh nodes were to be moved off the image domain." );
	}

	this->m_DisplacementsUpdated = true;
	this->m_RegionsUpdated = (changed==0);
	this->m_EnergyUpdated = (changed==0);
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	if ( !this->m_EnergyUpdated ) {
		this->m_Value = 0.0;

		double vxvol = 1.0;
		for(size_t i = 0; i<Dimension; i++)
			vxvol *= this->GetCurrentMap(0)->GetSpacing()[i];

		size_t nPix = this->GetCurrentMap(0)->GetLargestPossibleRegion().GetNumberOfPixels();
		size_t nrois = m_ROIs.size();

		const ReferencePixelType* refBuffer = this->m_ReferenceImage->GetBufferPointer();
		const typename ProbabilityMapType::PixelType* bgBuffer = this->m_BackgroundMask->GetBufferPointer();
		const typename ProbabilityMapType::PixelType* roiBuffer[nrois];

		for( size_t roi = 0; roi < nrois; roi++ ) {
			roiBuffer[roi] = this->GetCurrentMap(roi)->GetBufferPointer();
		}

		ReferencePointType pos;
		ReferencePixelType val;
		typename ProbabilityMapType::PixelType w;
		typename ProbabilityMapType::PixelType bgw;
		double totalVol;
		MeasureType smpl_val;
		for( size_t i = 0; i < nPix; i++) {
			totalVol = 0.0;
			smpl_val = 0.0;
			bgw = *(bgBuffer + i);

			if (bgw >= 0.9) {
				continue;
			}

			val = *(refBuffer+i);
			for( size_t roi = 0; roi < nrois; roi++ ) {
				w = *( roiBuffer[roi] + i );
				if ( w > 1.0e-8 ) {
					smpl_val +=  w * this->GetEnergyOfSample( val, roi, true );
					totalVol += w;
				}
			}

			if (bgw > 1.0e-3) {
				smpl_val+= bgw * this->m_MaxEnergy;
				totalVol+= bgw;
			}

			this->m_Value+= smpl_val / totalVol;
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


// Reorient contours to image direction in order that allowing pixel-wise computations
// ReorientFilter computes the new extent of the image if the directions
// matrix is identity. This is necessary to be able to binarize the contours
// (that are given in physical coordinates).
// See https://github.com/oesteban/ACWE-Registration/issues/92
template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeCurrentContours() {
	// Copy contours
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		ContourCopyPointer copy = ContourCopyType::New();
		copy->SetInput( this->m_Priors[contid] );
		copy->Update();
		this->m_CurrentContours.push_back( copy->GetOutput() );

		ShapeCopyPointer copyShape = ShapeCopyType::New();
		copyShape->SetInput( this->m_Priors[contid] );
		copyShape->Update();
		this->m_Gradients.push_back( copyShape->GetOutput() );
	}

	// Fill in interpolator points
	this->m_CurrentDisplacements = PointDataContainer::New();
	this->m_CurrentDisplacements->Reserve( this->m_NumberOfPoints );
	VectorType zerov;
	zerov.Fill(0.0);

	size_t identifier = 0;
	for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
		PointsIterator c_it  = this->m_CurrentContours[contid]->GetPoints()->Begin();
		PointsIterator c_end = this->m_CurrentContours[contid]->GetPoints()->End();
		PointType ci;
		while( c_it != c_end ) {
			ci = c_it.Value();
			this->m_NodesPosition.push_back( ci );
			this->m_CurrentDisplacements->SetElement(identifier, zerov );
			++c_it;
			identifier++;
		}
	}

	if ( this->m_NumberOfPoints != this->m_NodesPosition.size() ) {
		itkExceptionMacro( << "an error occurred initializing mesh points: NumberOfPoints in functional and" \
				" NumberOfSamples in transform do not match" );

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
::ComputeOuterRegions() {
	ContourOuterRegions outerVect;

	if( this->m_NumberOfRegions > 2 ) {
		// Set up ROI interpolator
		typename ROIInterpolatorType::Pointer interp = ROIInterpolatorType::New();
		interp->SetInputImage( this->m_CurrentRegions );

		// Set up outer regions
		for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
			// Compute mesh of normals
			NormalFilterPointer normalsFilter = NormalFilterType::New();
			normalsFilter->SetInput( this->m_CurrentContours[contid] );
			normalsFilter->Update();
			ContourPointer normals = normalsFilter->GetOutput();
			outerVect.resize( normals->GetNumberOfPoints() );

			typename ContourType::PointsContainerConstIterator c_it  = normals->GetPoints()->Begin();
			typename ContourType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

			PointType ci;
			VectorType v;
			VectorType ni;

			size_t pid;
			while( c_it != c_end ) {
				pid = c_it.Index();
				ci = c_it.Value();
				normals->GetPointData( pid, &ni );
				ROIPixelType inner = interp->Evaluate( ci + ni );
				ROIPixelType outer = interp->Evaluate( ci - ni );
				outerVect[pid] = outer;
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
double
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputePointArea(const PointIdentifier & iId, ContourType *mesh ) {
	QEType* edge = mesh->FindEdge( iId );
	QEType* temp = edge;
	CellIdentifier cell_id(0);
	double totalArea = 0.0;
	ContourPointType pt[3];
	typedef typename PolygonType::PointIdIterator PolygonPointIterator;

	do {
		cell_id = temp->GetLeft();

		if ( cell_id != ContourType::m_NoFace ) {
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

	if ( npoints != this->m_NumberOfPoints ) {
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
		this->m_OldDirection = _arg->GetDirection();
		this->m_OldOrigin = _arg->GetOrigin();

		DirectionType itk; itk.SetIdentity();
		itk(0,0) = -1.0; itk(1,1) = -1.0;
		if ( this->m_OldDirection != itk ) {
			itkWarningMacro( << "input volume is not RAS oriented.")
		}

		//typedef itk::OrientImageFilter< ReferenceImageType, ReferenceImageType >  Orienter;
		//typename Orienter::Pointer orient = Orienter::New();
		//orient->UseImageDirectionOn();
		//orient->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
		//orient->SetInput(_arg);
		//orient->Update();
		//ReferenceImagePointer ref = orient->GetOutput();

		ReferenceImagePointer ref = ReferenceImageType::New();
		ref->CopyInformation( _arg );
		ref->SetRegions( _arg->GetLargestPossibleRegion() );
		ref->Allocate();

		itk::ImageAlgorithm::Copy< ReferenceImageType, ReferenceImageType >(
				_arg,
				ref,
				_arg->GetLargestPossibleRegion(),
				ref->GetLargestPossibleRegion()
		);


		DirectionType dir; dir.SetIdentity();
		ref->SetDirection( dir );
		ref->SetOrigin( PointType(0.0) + this->m_OldDirection*(this->m_OldOrigin).GetVectorFromOrigin() );

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
::EvaluateGradient(  typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point,
		size_t outer_roi, size_t inner_roi ) const {
	if ( outer_roi == inner_roi ) {
		return 0.0;
	}
	ReferencePixelType value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	MeasureType grad = this->GetEnergyOfSample( value, outer_roi ) - this->GetEnergyOfSample( value, inner_roi ) + isOutside * this->m_MaxEnergy;
	grad = (fabs(grad)>MIN_GRADIENT)?grad:0.0;
	return grad;
}


template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetEnergyAtPoint( typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point, size_t roi ) const {
	ReferencePixelType value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	return this->GetEnergyOfSample( value, roi ) + isOutside * this->m_MaxEnergy;
}

template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetEnergyAtPoint( typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point, size_t roi,
		            typename FunctionalBase<TReferenceImageType, TCoordRepType>::ReferencePixelType & value) const {
	value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	return this->GetEnergyOfSample( value, roi ) + isOutside * this->m_MaxEnergy;
}

}

#endif /* FUNCTIONALBASE_HXX_ */
