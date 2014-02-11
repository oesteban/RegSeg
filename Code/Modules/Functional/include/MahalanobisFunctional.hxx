// --------------------------------------------------------------------------
// File:             MahalanobisFunctional.hxx
// Date:             27/10/2012
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
// This file is part of ACWERegistration-Debug@Debug
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

#ifndef MAHALANOBISLEVELSETS_HXX_
#define MAHALANOBISLEVELSETS_HXX_

#include "MahalanobisFunctional.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_ldl_cholesky.h>

#include "DisplacementFieldFileWriter.h"
#include <itkMeshFileWriter.h>

namespace rstk {

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
::PrintSelf( std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "NumberOfRegions: " << this->m_NumberOfRegions << std::endl;
	for( size_t roi = 0; roi < this->m_NumberOfRegions; roi++ ) {
		os << indent << "Region " << roi << ":" << std::endl;
		os << indent << indent << "Mean vector = " << this->m_Parameters[roi].mean << std::endl;
		os << indent << indent << "Covariance Matrix = " << std::endl << this->m_Parameters[roi].cov << std::endl;
		os << indent << indent << "Covariance Matrix^-1 = " << std::endl << this->m_Parameters[roi].invcov << std::endl;
	}
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
::Initialize() {
	Superclass::Initialize();

	// Check that parameters are initialized
	if (! this->ParametersInitialized() )
		this->UpdateDescriptors();
}

template <typename TReferenceImageType, typename TCoordRepType>
inline typename MahalanobisFunctional<TReferenceImageType,TCoordRepType>::MeasureType
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
::GetEnergyOfSample( typename MahalanobisFunctional<TReferenceImageType,TCoordRepType>::ReferencePixelType value, size_t roi ) const {
	ReferencePixelType dist = value - this->m_Parameters[roi].mean;
	return dot_product(dist.GetVnlVector(), this->m_Parameters[roi].invcov.GetVnlMatrix() * dist.GetVnlVector() ) + this->m_Parameters[roi].bias;
}


template <typename TReferenceImageType, typename TCoordRepType>
void MahalanobisFunctional<TReferenceImageType,TCoordRepType>
::UpdateDescriptors() {
	// Update regions
	for( size_t roi = 0; roi < this->m_NumberOfRegions; roi++ ) {
		ParametersType param = this->UpdateParametersOfRegion(roi);
		this->SetParameters(roi, param);
	}
}

template <typename TReferenceImageType, typename TCoordRepType>
typename MahalanobisFunctional<TReferenceImageType,TCoordRepType>::ParametersType
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
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
	typename CovarianceFilter::MatrixType cov = covFilter->GetCovarianceMatrix();

	for (size_t row = 0; row < Components; row ++ ) {
		for (size_t col = 0; col < Components; col ++ ) {
			newParameters.cov(row,col) = cov(row,col);
		}
	}

	return newParameters;
}




template <typename TReferenceImageType, typename TCoordRepType>
size_t
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
::AddShapePrior( typename MahalanobisFunctional<TReferenceImageType,TCoordRepType>::ContourType* prior,
		         typename MahalanobisFunctional<TReferenceImageType,TCoordRepType>::ParametersType& params){
	size_t id = this->AddShapePrior( prior );
	this->SetParameters( id, params );
	return id;
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
::SetParameters( size_t roi,
		         typename MahalanobisFunctional<TReferenceImageType,TCoordRepType>::ParametersType& params ) {

	this->m_Parameters[roi].mean = params.mean;

	CovarianceType cov = params.cov;

	if (cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols()) {
		itkExceptionMacro(<< "Covariance matrix must be square");
	}

	if ( cov.GetVnlMatrix().rows() != Components ) {
		itkExceptionMacro(<< "Length of measurement vectors must be the same as the size of the covariance.");
	}

	ReferenceValueType det = 0.0;

	if( Components > 1 ) {
		// Compute diagonal and check that eigenvectors >= 0.0
		typedef typename vnl_diag_matrix<ReferenceValueType>::iterator DiagonalIterator;
		typedef vnl_symmetric_eigensystem<ReferenceValueType> Eigensystem;
		vnl_matrix< ReferenceValueType > vnlCov = cov.GetVnlMatrix();
		Eigensystem* e = new Eigensystem( vnlCov );

		bool modified = false;
		DiagonalIterator itD = e->D.begin();
		while ( itD!= e->D.end() ) {
			if (*itD < 0) {
				*itD = 0.;
				modified = true;
			}
			itD++;
		}

		if (modified)
			this->m_Parameters[roi].cov = e->recompose();
		else
			this->m_Parameters[roi].cov = cov;
		delete e;

		// the inverse of the covariance matrix is first computed by SVD
		vnl_matrix_inverse< ReferenceValueType > inv_cov( cov.GetVnlMatrix() );

		// the determinant is then costless this way
		det = inv_cov.determinant_magnitude();

		if( det < 0.) {
			itkExceptionMacro( << "| sigma | < 0" );
		}

		// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
		const ReferenceValueType singularThreshold = 1.0e-10;
		if( det > singularThreshold ) {
			// allocate the memory for inverse covariance matrix
			this->m_Parameters[roi].invcov = inv_cov.inverse();
		} else {
			// TODO Perform cholesky diagonalization and select the semi-positive aproximation
			vnl_matrix< double > diag_cov( Components, Components );
			for ( size_t i = 0; i<Components; i++)
				for ( size_t j = 0; j<Components; j++)
					diag_cov[i][j] = vnlCov[i][j];
			vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( diag_cov );
			vnl_vector< double > D( chol->diagonal() );
			det = dot_product( D, D );
			vnl_matrix_inverse< double > R (chol->upper_triangle());

			for ( size_t i = 0; i<Components; i++)
				for ( size_t j = 0; j<Components; j++)
					this->m_Parameters[roi].invcov(i,j) = R.inverse()[i][j];
		}
	} else {
		this->m_Parameters[roi].cov = cov;
		this->m_Parameters[roi].invcov(0,0)=1.0 / cov(0,0);
		det = fabs( this->m_Parameters[roi].cov(0,0) );
	}
	this->m_Parameters[roi].bias = log( det );
	this->m_Parameters[roi].initialized = true;
}

template< typename TReferenceImageType, typename TCoordRepType >
bool
MahalanobisFunctional<TReferenceImageType, TCoordRepType>
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
MahalanobisFunctional<TReferenceImageType,TCoordRepType>
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
		ss << "], \"cov\": [ ";

		for( size_t j = 0; j<this->m_Parameters[i].cov.GetVnlMatrix().rows(); j++ ) {
			for( size_t k = 0; k<this->m_Parameters[i].cov.GetVnlMatrix().cols(); k++ ) {
				if( j>0 || k>0 ) ss << ",";
				ss << this->m_Parameters[i].cov(j,k);
			}
		}
		ss << "] }";
	}
	ss << "] } }";

	return ss.str();
}

}

#endif /* MAHALANOBISLEVELSETS_HXX_ */
