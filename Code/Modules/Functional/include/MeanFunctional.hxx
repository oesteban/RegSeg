// --------------------------------------------------------------------------------------
// File:          MeanFunctional.hxx
// Date:          Feb 10, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef MEANFUNCTIONAL_HXX_
#define MEANFUNCTIONAL_HXX_

#include "MeanFunctional.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_ldl_cholesky.h>

#include "DisplacementFieldFileWriter.h"
#include <itkMeshFileWriter.h>

namespace rstk {

template <typename TReferenceImageType, typename TCoordRepType>
void
MeanFunctional<TReferenceImageType,TCoordRepType>
::PrintSelf( std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "NumberOfRegions: " << this->m_NumberOfRegions << std::endl;
	for( size_t roi = 0; roi < this->m_NumberOfRegions; roi++ ) {
		os << indent << "Region " << roi << ":" << std::endl;
		os << indent << indent << "Mean vector = " << this->m_Parameters[roi].mean << std::endl;
	}
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MeanFunctional<TReferenceImageType,TCoordRepType>
::Initialize() {
	Superclass::Initialize();

	// Check that parameters are initialized
	if (! this->ParametersInitialized() )
		this->UpdateDescriptors();
}

template <typename TReferenceImageType, typename TCoordRepType>
inline typename MeanFunctional<TReferenceImageType,TCoordRepType>::MeasureType
MeanFunctional<TReferenceImageType,TCoordRepType>
::GetEnergyOfSample( typename MeanFunctional<TReferenceImageType,TCoordRepType>::ReferencePixelType value, size_t roi ) const {
	ReferencePixelType dist = value - this->m_Parameters[roi].mean;
	MeasureType val = dist.GetNorm();
	return val * val;
}

template <typename TReferenceImageType, typename TCoordRepType>
void MeanFunctional<TReferenceImageType,TCoordRepType>
::UpdateDescriptors() {
	// Update regions
	for( size_t roi = 0; roi < this->m_NumberOfRegions; roi++ ) {
		ParametersType param = this->UpdateParametersOfRegion(roi);
		this->SetParameters(roi, param);
	}
}

template <typename TReferenceImageType, typename TCoordRepType>
typename MeanFunctional<TReferenceImageType,TCoordRepType>::ParametersType
MeanFunctional<TReferenceImageType,TCoordRepType>
::UpdateParametersOfRegion( const size_t idx ) {
	ParametersType newParameters;
	ProbabilityMapConstPointer roipm = this->GetCurrentMap( idx );

	// Apply weighted mean/covariance estimators from ITK
	typename CovarianceFilter::WeightArrayType weights;
	ReferenceSamplePointer sample = ReferenceSampleType::New();
	sample->SetImage( this->m_ReferenceImage );

	size_t sampleSize = this->m_ReferenceImage->GetLargestPossibleRegion().GetNumberOfPixels();
	weights.SetSize( sampleSize );
	weights.Fill( 0.0 );

	const typename ProbabilityMapType::PixelType* roipmb = roipm->GetBufferPointer();
	for ( size_t pidx = 0; pidx < sampleSize; pidx++) {
		weights[pidx] = *( roipmb + pidx );
	}

	CovarianceFilterPointer covFilter = CovarianceFilter::New();
	covFilter->SetInput( sample );
	covFilter->SetWeights( weights );
	covFilter->Update();

	newParameters.mean = covFilter->GetMean();
	return newParameters;
}




template <typename TReferenceImageType, typename TCoordRepType>
size_t
MeanFunctional<TReferenceImageType,TCoordRepType>
::AddShapePrior( typename MeanFunctional<TReferenceImageType,TCoordRepType>::ContourType* prior,
		         typename MeanFunctional<TReferenceImageType,TCoordRepType>::ParametersType& params){
	size_t id = this->AddShapePrior( prior );
	this->SetParameters( id, params );
	return id;
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MeanFunctional<TReferenceImageType,TCoordRepType>
::SetParameters( size_t roi,
		         typename MeanFunctional<TReferenceImageType,TCoordRepType>::ParametersType& params ) {

	this->m_Parameters[roi].mean = params.mean;
	this->m_Parameters[roi].initialized = true;
}

template< typename TReferenceImageType, typename TCoordRepType >
bool
MeanFunctional<TReferenceImageType, TCoordRepType>
::ParametersInitialized() {
	if ( this->m_Parameters.size() != this->m_NumberOfRegions ) {
		this->m_Parameters.resize( this->m_NumberOfRegions );
		for ( size_t i = 0; i<this->m_NumberOfRegions; i++ )
			this->m_Parameters[i].initialized = false;
		return false;
	}

	for ( size_t i = 0; i < this->m_NumberOfRegions; i++ ) {
		if ( ! this->m_Parameters[i].initialized ) return false;
	}
	return true;
}

template <typename TReferenceImageType, typename TCoordRepType>
std::string
MeanFunctional<TReferenceImageType,TCoordRepType>
::PrintFormattedDescriptors() {
	std::stringstream ss;

	ss << "{ \"descriptors\" : { \"number\": " << this->m_NumberOfRegions << ", \"values\": [";

	for ( size_t i = 0; i<this->m_NumberOfRegions; i++ ){
		if (i>0) ss<<",";

		ss << "{ \"id\": " << i << ", \"mu\": [";

		for ( size_t l = 0; l<this->m_Parameters[i].mean.Size(); l++ ) {
			if( l>0 ) ss << ",";
			ss << this->m_Parameters[i].mean[l];
		}

		ss << "] }";
	}
	ss << "] } }";

	return ss.str();
}

}
#endif /* MEANFUNCTIONAL_HXX_ */
