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

#include "DisplacementFieldFileWriter.h"

namespace rstk {
template <typename TReferenceImageType, typename TCoordRepType>
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::MahalanobisLevelSets() {


}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::PrintSelf( std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	// TODO PrintSelf parameters
}

template <typename TReferenceImageType, typename TCoordRepType>
typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::MeasureType
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::GetValue() const {
	// for all classes

		// if pixel is inside object

			// compute mahalanobis distance to class

			// add to energy

	// reset m_Value and return
	return this->m_Value;
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::SetParameters( typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::MeanType& mean,
		         typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::CovarianceType& cov,
		         bool inside ) {
	unsigned int idx = (unsigned int ) (!inside);

	this->m_Mean[idx] = mean;

	if (cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols()) {
		itkExceptionMacro(<< "Covariance matrix must be square");
	}

	if ( cov.GetVnlMatrix().rows() != Components ) {
		itkExceptionMacro(<< "Length of measurement vectors must be the same as the size of the covariance.");
	}

	if( Components > 1 ) {
		// Compute diagonal and check that eigenvectors >= 0.0
		typedef typename vnl_diag_matrix<PixelValueType>::iterator DiagonalIterator;
		typedef vnl_symmetric_eigensystem<PixelValueType> Eigensystem;
		vnl_matrix< PixelValueType > S = cov.GetVnlMatrix();
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
		vnl_matrix_inverse< PixelValueType > inv_cov( this->m_InverseCovariance[idx].GetVnlMatrix() );

		// the determinant is then costless this way
		PixelValueType det = inv_cov.determinant_magnitude();

		if( det < 0.) {
			itkExceptionMacro( << "det( m_InverseCovariance ) < 0" );
		}

		// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
		const PixelValueType singularThreshold = 1.0e-6;
		if( det > singularThreshold ) {
			// allocate the memory for m_InverseCovariance matrix
			this->m_InverseCovariance[idx].GetVnlMatrix() = inv_cov.inverse();
		} else {
			// TODO Perform cholesky diagonalization and select the semi-positive aproximation
			//vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( this->m_InverseCovariance[idx].GetVnlMatrix() );
			//vnl_vector< PixelValueType > D( chol->diagonal() );
			//det = dot_product( D, D );
			//vnl_matrix_inverse< PixelValueType > R (chol->upper_triangle());
			//this->m_InverseCovariance[idx].GetVnlMatrix() = R.inverse();
		}
	} else {
		this->m_InverseCovariance[idx](0,0)=1.0/cov(0,0);
	}
}

template <typename TReferenceImageType, typename TCoordRepType>
typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::DeformationFieldPointer
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::GetLevelSetsMap( MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::DeformationFieldType* levelSetMap) {
	// Initialize interpolator
	InterpolatorPointer interp = InterpolatorType::New();
	interp->SetInputImage( this->m_ReferenceImage );

	// Compute mesh of normals
	NormalFilterPointer normFilter = NormalFilterType::New();
	normFilter->SetInput( this->m_CurrentContourPosition );
	normFilter->Update();
	ContourDeformationPointer normals = normFilter->GetOutput();

	typename ContourDeformationType::PointsContainerConstIterator c_it = normals->GetPoints()->Begin();
	typename ContourDeformationType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

	PointValueType sign[2] = { -1.0, 1.0 };  // IMPORTANT: signs are reversed, normal should be inwards but filter
	                                         // gives outwards.
	PointValueType levelSet;
	PointType  ci;          //
	PointType  ci_prime;
	PixelType  fi;          // Feature on ci_prime
	typename ContourDeformationType::PixelType ni;
	typename ContourDeformationType::PointIdentifier idx;

	// TODO: for all existing contours

	// for all node in mesh
	while (c_it!=c_end) {
		levelSet = 0;
		idx = c_it.Index();
		ci_prime = c_it.Value();
		fi = interp->Evaluate( ci_prime );    // Feature in c'_{i}
		for( size_t i = 0; i<2; i++) {                 // Compute on both sides of the levelset
			PixelType dist = fi - m_Mean[i];
			// compute mahalanobis distance in position
			levelSet+= sign[i] * dot_product(dist.GetVnlVector(), m_InverseCovariance[i].GetVnlMatrix() * dist.GetVnlVector() );
		}
		assert( !std::isnan(levelSet) );
	    // project to normal, updating transform
		normals->GetPointData( idx, &ni );         // Normal ni in point c'_i
		ni*=levelSet;
		normals->SetPointData( idx, ni );
		++c_it;
	}

	// Interpolate sparse velocity field to targetDeformation
	this->m_SparseToDenseResampler->CopyImageInformation( levelSetMap );
	this->m_SparseToDenseResampler->SetInput( normals );
	this->m_SparseToDenseResampler->Update();
	return this->m_SparseToDenseResampler->GetOutput();

	/*
	// Output map (testing purposes)
	typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;
	typename Writer::Pointer writer = Writer::New();
	writer->SetFileName( "speedtest.nii.gz" );
	writer->SetInput( levelSetMap );
	writer->Update();*/

	/*
	// Write final result out
	typedef itk::VTKPolyDataWriter< ContourDeformationType >     WriterType;
	typename WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( normals );
	polyDataWriter->SetFileName( "speedmap.vtk" );
	polyDataWriter->Update();
	*/

}

}

#endif /* MAHALANOBISLEVELSETS_HXX_ */
