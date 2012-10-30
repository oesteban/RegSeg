// --------------------------------------------------------------------------
// File:             MahalanobisLevelSets.hxx
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

#include "MahalanobisLevelSets.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_ldl_cholesky.h>

namespace rstk {
template <class TTargetImage, class TDeformationField, class TContourDeformation>
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::MahalanobisLevelSets() {


}

template <class TTargetImage, class TDeformationField, class TContourDeformation>
typename MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::ValueType
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::GetValue() const {
	// for all classes

		// if pixel is inside object

			// compute mahalanobis distance to class

			// add to energy

	// reset m_Value and return
	return this->m_Value;
}

template <class TTargetImage, class TDeformationField, class TContourDeformation>
void
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::SetParameters( typename MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::MeanType& mean,
		         typename MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::CovarianceType& cov,
		         bool inside ) {
	unsigned int idx = (unsigned int ) (!inside);

	this->m_Mean[idx] = mean;

	// TODO Check covariance size
	if (cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols()) {
		itkExceptionMacro(<< "Covariance matrix must be square");
	}

/*	if ( this->GetMeasurementVectorSize() ){
		if ( cov.GetVnlMatrix().rows() != this->GetMeasurementVectorSize() ) {
		itkExceptionMacro(<< "Length of measurement vectors must be the same as the size of the covariance.");
		}
	} else {
		// not already set, cache the size
		this->SetMeasurementVectorSize( cov.GetVnlMatrix().rows() );
	}

	// TODO
	if (m_InverseCovariance == cov) { // no need to copy the matrix, compute the inverse, or the normalization
		return;
	}
*/

	// Compute diagonal and check that eigenvectors >= 0.0
	typedef vnl_diag_matrix<double>::iterator DiagonalIterator;
	typedef vnl_symmetric_eigensystem<double> Eigensystem;
	vnl_matrix< double > S = cov.GetVnlMatrix();
	Eigensystem* e = new Eigensystem( S );

	bool modified = false;
	DiagonalIterator itD = e->D.begin();
	while ( itD!= e->D.end() ) {
		if (*itD < 0) {
			*itD = 0.;
			modified = true;
		}
		itD++;
	}

	if (modified) this->m_InverseCovariance[idx] = e->recompose();
	else this->m_InverseCovariance[idx] = cov;
	delete e;

	// the inverse of the covariance matrix is first computed by SVD
	vnl_matrix_inverse< double > inv_cov( this->m_InverseCovariance[idx].GetVnlMatrix() );

	// the determinant is then costless this way
	double det = inv_cov.determinant_magnitude();

	if( det < 0.) {
		itkExceptionMacro( << "det( m_InverseCovariance ) < 0" );
	}

	// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
	const double singularThreshold = 1.0e-6;
	if( det > singularThreshold ) {
		// allocate the memory for m_InverseCovariance matrix
		this->m_InverseCovariance[idx].GetVnlMatrix() = inv_cov.inverse();
	} else {
		// Perform cholesky diagonalization and select the semi-positive aproximation
		vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( this->m_InverseCovariance[idx].GetVnlMatrix() );
		vnl_vector< double > D( chol->diagonal() );
		det = dot_product( D, D );
		vnl_matrix_inverse< double > R (chol->upper_triangle());
		this->m_InverseCovariance[idx].GetVnlMatrix() = R.inverse();
	}
}

template <class TTargetImage, class TDeformationField, class TContourDeformation>
void
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::GetLevelSetsMap( MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::DeformationFieldType & levelSetMap) const {
	// Copy deformation map
	ContourDeformationPointer speedMap = this->m_ContourDeformation->Copy();

	// Compute mesh of normals
	NormalFilterPointer normFilter = NormalFilterType::New();
	normFilter->SetInput( this->m_ContourDeformation );
	normFilter->Update();
	ContourDeformationPointer normals = normFilter->GetOutput();

	typename ContourDeformationType::PointDataContainerPointer points = this->m_ContourDeformation->GetPoints();
	typename ContourDeformationType::PointDataContainerIterator p_it = points->Begin();
	typename ContourDeformationType::PointDataContainerPointer container = this->m_ContourDeformation->GetPointData();
	typename ContourDeformationType::PointDataContainerIterator u_it = container->Begin();
	typename ContourDeformationType::PointDataContainerPointer normContainer = normals->GetPointData();
	typename ContourDeformationType::PointDataContainerIterator n_it = normContainer->Begin();
	//typename ContourDeformationType::PointDataContainerPointer speeds = speedMap->GetPointData();
	//typename ContourDeformationType::PointDataContainerIterator s_it = speeds->Begin();

	// for all node in mesh
	while (p_it != container->End()) {
		VectorValueType levelSet = 0;
		ValueType uk = u_it.Value();
		ValueType currentPoint = p_it.Value() + uk;
		ValueType Ik = this->m_Image->GetValue( currentPoint );
		// compute on both segments
		for( size_t i = 0; i<2; i++) {
			ValueType Dk = Ik - m_Mean[i];
			// compute mahalanobis distance in position
			levelSet-= dot_product(Dk, m_InverseCovariance[i].GetVnlMatrix() * Dk);
		}
	    // project to normal, updating transform
		speedMap->SetPointData( p_it.Index(), levelSet * n_it.Value() );
		++p_it;
		++u_it;
		++n_it;
	}

	// Interpolate sparse velocity field to targetDeformation
	typename ResamplerType::Pointer res = ResamplerType::New();
	res->SetReferenceImage( this->m_DeformationField );
	res->SetInput( speedMap );
	res->Update();

	levelSetMap = res->GetOutput()->Copy();
}

}

#endif /* MAHALANOBISLEVELSETS_HXX_ */
