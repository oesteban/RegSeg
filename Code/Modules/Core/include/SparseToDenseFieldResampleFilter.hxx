// --------------------------------------------------------------------------
// File:             SparseToDenseFieldResampleFilter.hxx
// Date:             30/10/2012
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

#ifndef SPARSETODENSEFIELDRESAMPLEFILTER_HXX_
#define SPARSETODENSEFIELDRESAMPLEFILTER_HXX_

#include "SparseToDenseFieldResampleFilter.h"
#include "IDWMultivariateInterpolator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace rstk {

template<class TInputMesh, class TOutputImage>
SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::SparseToDenseFieldResampleFilter() {
	m_OutputSpacing.Fill(1.0);
	m_OutputOrigin.Fill(0.0);
	m_OutputDirection.SetIdentity();
	m_OutputSize.Fill(0);
	m_OutputStartIndex.Fill(0);

	m_Interpolator =
			IDWMultivariateInterpolator<InputMeshType, OutputImageType>::New();
	m_DefaultPixelValue = NumericTraits<ValueType>::ZeroValue(
			m_DefaultPixelValue);
	m_IsPhiInitialized = false;
	m_N = 0;
	m_k = 0;
};

template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::
AddControlPoints( const typename SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::InputMeshType* prior) {
	for( size_t i = 0; i<prior->GetNumberOfPoints();i++ ) {
		this->m_ControlPoints.push_back( prior->GetPoint(i));
	}
	this->m_N = this->m_ControlPoints.size();
}

template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::PrintSelf(
		std::ostream &os, Indent indent) const {
	//Superclass::PrintSelf(os, indent);

	os << indent << "DefaultPixelValue: "
			<< static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_DefaultPixelValue)
			<< std::endl;
	os << indent << "Size: " << m_OutputSize << std::endl;
	os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
	os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
	os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
	os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
	os << indent << "Interpolator: " << m_Interpolator.GetPointer()
			<< std::endl;
}
;

/**
 * Set the output image spacing.
 */
template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::SetOutputSpacing(
		const double *spacing) {
	OutputSpacingType s(spacing);

	this->SetOutputSpacing(s);
}

/**
 * Set the output image origin.
 */
template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::SetOutputOrigin(
		const double *origin) {
	OutputPointType p(origin);
	this->SetOutputOrigin(p);
}

/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::GenerateData() {
	if ( m_N == 0 ) {
		itkExceptionMacro(<< "Shape priors not set");
	}

	if (m_DefaultPixelValue.Size() == 0) {
		NumericTraits<OutputPixelType>::SetLength(m_DefaultPixelValue,
				this->GetOutput()->GetNumberOfComponentsPerPixel());
		m_DefaultPixelValue.Fill(0);
	}

//	this->GetOutput()->FillBuffer( this->m_DefaultPixelValue );

	m_k = this->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();

	for (size_t i = 0; i < OutputImageDimension; i++ ) m_LevelSetVector[i] = SpeedsVector( m_N );
	size_t nConts = this->GetNumberOfIndexedInputs();

	if (!m_IsPhiInitialized) {
		m_Phi = WeightsMatrix(m_k, m_N);
		double dist, wi;
		OutputPointType ci, pi;
		size_t row, col;

		// Walk the output region
		for( row = 0; row < this->m_k; row++ ) {
			this->GetOutput()->TransformIndexToPhysicalPoint(this->GetOutput()->ComputeIndex(row), pi);

			for ( col=0; col < this->m_N; col++) {
				// TODO Extract RBF from here (use a functor?)
				// Compute weight i: wi(x) (inverse distance)
				dist = (this->m_ControlPoints[col] - pi).GetNorm();
				wi = (dist < vnl_math::eps) ? 1.0 : (1.0 / vcl_pow(dist, 2));
				if (wi > 1e-3) {
					m_Phi.put(row, col, wi);
				}
			}
		}

		m_IsPhiInitialized = true;
	}

	size_t contsOffset = 0;
	for( size_t cont = 0; cont< nConts; cont++ ) {
		OutputPixelType ni;
		size_t nPoints = this->GetInput(cont)->GetNumberOfPoints();
		for (size_t el = 0; el < nPoints; el++) {
			this->GetInput(cont)->GetPointData( el, &ni );

			for ( size_t i = 0; i < OutputImageDimension; i++ ) {
				if( std::isnan( ni[i] )) ni[i] = 0;
				m_LevelSetVector[i].put( contsOffset + el, ni[i] );
			}
		}
		contsOffset += nPoints;
	}

	double norms = 0.0;
	for (size_t i = 0; i < OutputImageDimension; i++ ) {
		norms+= m_LevelSetVector[i].one_norm();
	}


	if( norms==0 ) return;

    SpeedsVector outvector[3];
    norms=0;
    for ( size_t i = 0; i < OutputImageDimension; i++ ) {
    	m_Phi.mult(m_LevelSetVector[i], outvector[i] );
    	norms+= outvector[i].one_norm();
    }

    if( norms==0 ) return;

    bool isZero;
    OutputPixelType* buffer = this->GetOutput()->GetBufferPointer();
    OutputPixelType ni;
    for( size_t k = 0; k<m_k; k++ ) {
    	ni = m_DefaultPixelValue;
    	isZero = true;
    	for ( size_t i = 0; i < OutputImageDimension; i++ ) {
    		ni[i] = outvector[i].get(k);
    		isZero = isZero && (ni[i]==0.0);
    	}
    	if ( isZero ) continue;
    	*( buffer + k ) = ni;
    }

}

template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::CopyImageInformation(
		const OutputImageType* image) {
	m_OutputDirection = image->GetDirection();
	m_OutputOrigin = image->GetOrigin();
	m_OutputSpacing = image->GetSpacing();
	m_OutputStartIndex = image->GetRequestedRegion().GetIndex();
	m_OutputSize = image->GetRequestedRegion().GetSize();
}

/**
 * Inform pipeline of required output region
 */
template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::UpdateOutputInformation() {
	this->Modified();
	this->Superclass::UpdateOutputInformation();
}

/**
 * Inform pipeline of required output region
 */
template<class TInputMesh, class TOutputImage>
void SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::GenerateOutputInformation() {
	// call the superclass' implementation of this method
	this->Superclass::GenerateOutputInformation();

	// get pointers to the input and output
	OutputImagePointer outputPtr = this->GetOutput();
	if (!outputPtr) {
		return;
	}

	// Set the size of the output region
	typename TOutputImage::RegionType outputLargestPossibleRegion;
	outputLargestPossibleRegion.SetSize(m_OutputSize);
	outputLargestPossibleRegion.SetIndex(m_OutputStartIndex);
	outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
	outputPtr->SetBufferedRegion(outputLargestPossibleRegion);
	outputPtr->Allocate();
	outputPtr->FillBuffer(m_DefaultPixelValue);

	// Set spacing and origin
	outputPtr->SetSpacing(m_OutputSpacing);
	outputPtr->SetOrigin(m_OutputOrigin);
	outputPtr->SetDirection(m_OutputDirection);

	return;
}

/**
 * Verify if any of the components has been modified.
 */
template<class TInputMesh, class TOutputImage>
unsigned long SparseToDenseFieldResampleFilter<TInputMesh, TOutputImage>::GetMTime(
		void) const {
	unsigned long latestTime = Object::GetMTime();
/*
	if (this->GetInput().IsNotNull()) {
		if (latestTime < this->GetInput()->GetMTime()) {
			latestTime = this->GetInput()->GetMTime();
		}
	}*/
	return latestTime;
}

}

#endif /* SPARSETODENSEFIELDRESAMPLEFILTER_HXX_ */
