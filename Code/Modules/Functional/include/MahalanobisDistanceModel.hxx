// --------------------------------------------------------------------------------------
// File:          MahalanobisDistanceModel.hxx
// Date:          Dec 23, 2014
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
// This file is part of ACWEReg
//
// ACWEReg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWEReg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWEReg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef _MAHALANOBISDISTANCEMODEL_HXX_
#define _MAHALANOBISDISTANCEMODEL_HXX_

#include "MahalanobisDistanceModel.h"


#include <math.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_ldl_cholesky.h>


namespace rstk {
template< typename TInputVectorImage, typename TPriorsPrecisionType >
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::MahalanobisDistanceModel():
  Superclass(), m_MaxEnergyGap(0.0) {
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::PrintSelf(std::ostream & os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::InitializeMemberships() {
	this->m_NumberOfRegions = this->GetPriorsMap()->GetNumberOfComponentsPerPixel();
	this->m_Memberships.resize(this->m_NumberOfRegions);
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::GenerateData() {
	this->InitializeMemberships();

	ReferenceSamplePointer sample = ReferenceSampleType::New();
	sample->SetImage( this->GetInput() );
	size_t npix = sample->Size();
	size_t nregions = this->m_NumberOfRegions;

	std::vector<WeightArrayType> weights;
	const PriorsPrecisionType* priors = this->GetPriorsMap()->GetBufferPointer();

	for( size_t roi = 0; roi < nregions; roi++ ) {
		WeightArrayType w;
		w.SetSize(npix);

		for( size_t i = 0; i < npix; i++ ) {
			w[i] = *(priors + nregions * i + roi);
		}
		weights.push_back(w);
	}

	this->m_Means.empty();
	this->m_Covariances.empty();
	this->m_RegionOffsetContainer.SetSize(nregions);
	this->m_RegionOffsetContainer.Fill(0.0);

	for( size_t roi = 0; roi < nregions; roi++ ) {
		CovarianceFilterPointer covFilter = CovarianceFilter::New();
		covFilter->SetInput( sample );
		covFilter->SetWeights( weights[roi] );
		covFilter->Update();

		InternalFunctionPointer mf = InternalFunctionType::New();
		mf->SetMean( covFilter->GetMean() );

		CovarianceMatrixType cov = covFilter->GetCovarianceMatrix();
		mf->SetCovariance( cov );

		this->m_Memberships[roi] = mf;
		this->m_RegionOffsetContainer[roi] = log(2.0 * vnl_math::pi) * cov.Rows() + this->ComputeCovarianceDeterminant(cov);
		this->m_Means.push_back(covFilter->GetMean());
		this->m_Covariances.push_back(cov);
	}

	this->ComputeMaxEnergyGap();
}


template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::ComputeMaxEnergyGap() {
	MeasureType m;
	MeasureType v = itk::NumericTraits<MeasureType>::min();

	for( size_t roi1 = 0; roi1 < this->m_NumberOfRegions; roi1++) {
		for( size_t roi2 = 0; roi2 < this->m_NumberOfRegions; roi2++) {
			MeasurementVectorType mu = this->m_Means[roi1];
			m = this->Evaluate(mu, roi2);
			if (m > v)	v = m;
		}
	}
	this->m_MaxEnergyGap = fabs(v) * 1e6;
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
typename MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >::MeasureType
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::ComputeCovarianceDeterminant(CovarianceMatrixType& cov) const{
	typedef typename CovarianceMatrixType::ValueType ValueType;
	size_t ncomps = cov.Rows();
	CovarianceMatrixType invcov;
	MeasureType det = 1.0;

	if( ncomps > 1 ) {
		// Compute diagonal and check that eigenvectors >= 0.0
		typedef typename vnl_diag_matrix<ValueType>::iterator DiagonalIterator;
		typedef vnl_symmetric_eigensystem<ValueType> Eigensystem;
		vnl_matrix< ValueType > vnlCov = cov.GetVnlMatrix();
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
			cov = e->recompose();

		delete e;

		// the inverse of the covariance matrix is first computed by SVD
		vnl_matrix_inverse< ValueType > inv_cov( cov.GetVnlMatrix() );

		// the determinant is then costless this way
		det = inv_cov.determinant_magnitude();

		if( det < 0.) {
			itkExceptionMacro( << "| sigma | < 0" );
		}

		// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
		const ValueType singularThreshold = 1.0e-10;
		if( det > singularThreshold ) {
			// allocate the memory for inverse covariance matrix
			invcov = inv_cov.inverse();
		} else {
			// TODO Perform cholesky diagonalization and select the semi-positive aproximation
			vnl_matrix< double > diag_cov( ncomps, ncomps );
			for ( size_t i = 0; i<ncomps; i++)
				for ( size_t j = 0; j<ncomps; j++)
					diag_cov[i][j] = vnlCov[i][j];
			vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( diag_cov );
			vnl_vector< double > D( chol->diagonal() );
			det = dot_product( D, D );
			vnl_matrix_inverse< double > R (chol->upper_triangle());

			for ( size_t i = 0; i<ncomps; i++)
				for ( size_t j = 0; j<ncomps; j++)
					invcov(i,j) = R.inverse()[i][j];
		}
	} else {
		invcov(0,0)=1.0 / cov(0,0);
		det = fabs( cov(0,0) );
	}
	return log( det );
}

//template< typename TInputVectorImage, typename TPriorsPrecisionType >
//void
//MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
//::BeforeThreadedGenerateData() {
//	// find the actual number of threads
//	long nbOfThreads = this->GetNumberOfThreads();
//	if ( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 ) {
//		nbOfThreads = vnl_math_min( this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
//	}
//	// number of threads can be constrained by the region size, so call the
//	// SplitRequestedRegion
//	// to get the real number of threads which will be used
//	RegionType splitRegion;  // dummy region - just to call the following method
//	nbOfThreads = this->SplitRequestedRegion(0, nbOfThreads, splitRegion);
//
//	this->m_NumberOfRegions = this->GetPriorsMap()->GetNumberOfComponentsPerPixel();
//	this->InitializeMemberships();
//
//	this->m_ThreadMean.empty();
//	this->m_ThreadWeights.empty();
//	PixelType v;
//	v.Fill(0.0);
//
//	for( size_t t = 0; t < nbOfThreads; t++) {
//		MeanContainer mc;
//		WeightsContainer wc;
//		for( size_t i = 0; i < this->m_NumberOfRegions; i++ ) {
//			mc.push_back(v);
//			wc.push_back(0.0);
//		}
//		this->m_ThreadMean.push_back(mc);
//		this->m_ThreadWeights.push_back(wc);
//	}
//}
//
//template< typename TInputVectorImage, typename TPriorsPrecisionType >
//void
//MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
//::ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId) {
//	long nbOfPixels = inputRegionForThread.GetNumberOfPixels();
//	ProgressReporter progress( this, threadId, nbOfPixels );
//
//	MeanContainer mc = this->m_ThreadMean[threadId];
//	WeightsContainer wc = this->m_ThreadWeights[threadId];
//
//	ImageRegionConstIterator< TInputVectorImage > inputIt( this->GetInput(), inputRegionForThread );
//	ImageRegionConstIterator< PriorsImageType >   priorIt( this->GetPriorsMap(), inputRegionForThread );
//	ImageRegionConstIterator< MaskType >          maskIt ( this->GetMask(), inputRegionForThread );
// 	inputIt.GoToBegin();
//
// 	PriorsPixelType w;
// 	PriorsPrecisionType m;
// 	PixelType v;
// 	MeasureType e;
// 	v.Fill(0.0);
//
// 	while ( !inputIt.IsAtEnd() ) {
//		e = 0.0;
//		m = maskIt.Get();
//
//		if (m > 1.0e-5) {
//			v = inputIt.Get();
//			w = priorIt.Get();
//
//			for(size_t roi = 0; roi < this->m_NumberOfRegions; roi++ ) {
//				mc[roi]+= w[roi] * v;
//				wc[roi]+= w[roi];
//			}
//		}
//
//	    ++inputIt;
//	    ++priorIt;
//	    progress.CompletedPixel();
//	}
//}
//
//template< typename TInputVectorImage, typename TPriorsPrecisionType >
//void
//MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
//::AfterThreadedGenerateData() {
//
//	for(size_t roi = 0; roi < this->m_NumberOfRegions; roi++ ) {
//		PixelType m;
//		m.Fill(0.0);
//		PriorsPrecisionType w = 0.0;
//		for( size_t t = 0; t < this->m_ThreadMean.size(); t++) {
//			m+= this->m_ThreadMean[t][roi];
//			w+= this->m_ThreadWeights[t][roi];
//		}
//
//		InternalFunctionPointer mf = static_cast< InternalFunctionType* >(this->m_Memberships[roi]);
//		mf->SetMean(m/w);
//	}
//}


}


#endif /* _MAHALANOBISDISTANCEMODEL_HXX_ */
