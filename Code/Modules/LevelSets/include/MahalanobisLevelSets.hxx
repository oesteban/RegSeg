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
#include <itkMeshFileWriter.h>

namespace rstk {
template <typename TReferenceImageType, typename TCoordRepType>
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::MahalanobisLevelSets() {
	this->m_Interp = InterpolatorType::New();

}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::PrintSelf( std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	// TODO PrintSelf parameters
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::InitializeSamplingGrid() {
	this->m_ReferenceSamplingGrid = DeformationFieldType::New();
	typename ReferenceImageType::SpacingType sp = this->m_ReferenceImage->GetSpacing();
	double spacing = itk::NumericTraits< double >::max();

	for (size_t i = 0; i<Dimension; i++ ){
		if (sp[i] < spacing )
			spacing = sp[i];
	}

	sp.Fill( spacing * 0.25 );

	typename ReferenceImageType::PointType origin = this->m_ReferenceImage->GetOrigin();
	typename ReferenceImageType::SizeType size = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
	typename itk::ContinuousIndex<double, Dimension> idx;
	for (size_t i = 0; i<Dimension; i++ ){
		idx[i]=size[i]-1.0;
	}

	typename ReferenceImageType::PointType end;
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
	w->SetFileName( "test.nii.gz" );
	w->Update();
#endif
}

template <typename TReferenceImageType, typename TCoordRepType>
typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::MeasureType
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::GetValue() {
	this->m_Interp->SetInputImage( this->m_ReferenceImage );
	return Superclass::GetValue();
}

template <typename TReferenceImageType, typename TCoordRepType>
inline typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::MeasureType
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::GetEnergyAtPoint( typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::PixelPointType & point, size_t cont, size_t outside ) {
	PixelType dist = this->m_Interp->Evaluate( point ) - this->m_Parameters[cont].mean[outside];
	return dot_product(dist.GetVnlVector(), this->m_Parameters[cont].iCovariance[outside].GetVnlMatrix() * dist.GetVnlVector() );
}


template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::AddShapePrior( typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::ContourDeformationType* prior,
		         typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::ParametersType& params){
	this->Superclass::AddShapePrior( prior );

	size_t id = this->m_Parameters.size();
	this->m_Parameters.resize( id + 1 );
	this->SetParameters( id, params );

	if ( this->m_ReferenceImage.IsNotNull() ) {
		this->CheckExtents( prior );
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
bool
MahalanobisLevelSets<TReferenceImageType, TCoordRepType>
::CheckExtents( typename MahalanobisLevelSets<TReferenceImageType, TCoordRepType>::ContourDeformationType* prior ) const {
	typename ContourDeformationType::PointsContainerConstIterator u_it = prior->GetPoints()->Begin();
    typename ContourDeformationType::PointsContainerConstIterator u_end = prior->GetPoints()->End();

    PointType p;
    ContinuousIndex idx;
	while( u_it != u_end ) {
		if ( !this->m_ReferenceImage->TransformPhysicalPointToContinuousIndex(u_it.Value(), idx ) ) {
			p = u_it.Value();
			PointType origin;
			typename ReferenceImageType::IndexType tmp_idx;
			tmp_idx.Fill(0);
			this->m_ReferenceImage->TransformIndexToPhysicalPoint( tmp_idx, origin );
			PointType end;

			typename ReferenceImageType::SizeType size = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
			for ( size_t d=0; d< ReferenceImageType::ImageDimension; d++)
				tmp_idx[d] = size[d]  - 1;
			this->m_ReferenceImage->TransformIndexToPhysicalPoint( tmp_idx, end );
			itkExceptionMacro( << "Setting prior surface: vertex " << p << " is outside image extents (" << origin << ", " << end << ").");
			return false;
		}
		++u_it;
	}

	return true;
}

template <typename TReferenceImageType, typename TCoordRepType>
void
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::SetParameters( size_t contour_id,
		         typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::ParametersType& params ) {
	for ( size_t idx=0; idx<2; idx++ ) {
		this->m_Parameters[contour_id].mean[idx] = params.mean[idx];
		CovarianceType cov = params.iCovariance[idx];

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

			if (modified) this->m_Parameters[contour_id].iCovariance[idx] = e->recompose();
			else this->m_Parameters[contour_id].iCovariance[idx] = cov;
			delete e;

			// the inverse of the covariance matrix is first computed by SVD
			vnl_matrix_inverse< PixelValueType > inv_cov( this->m_Parameters[contour_id].iCovariance[idx].GetVnlMatrix() );

			// the determinant is then costless this way
			PixelValueType det = inv_cov.determinant_magnitude();

			if( det < 0.) {
				itkExceptionMacro( << "det( m_Parameters[contour_id].iCovariance ) < 0" );
			}

			// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
			const PixelValueType singularThreshold = 1.0e-6;
			if( det > singularThreshold ) {
				// allocate the memory for inverse covariance matrix
				this->m_Parameters[contour_id].iCovariance[idx].GetVnlMatrix() = inv_cov.inverse();
			} else {
				// TODO Perform cholesky diagonalization and select the semi-positive aproximation
				//vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( this->m_Parameters[contour_id].iCovariance[idx].GetVnlMatrix() );
				//vnl_vector< PixelValueType > D( chol->diagonal() );
				//det = dot_product( D, D );
				//vnl_matrix_inverse< PixelValueType > R (chol->upper_triangle());
				//this->m_Parameters[contour_id].iCovariance[idx].GetVnlMatrix() = R.inverse();
			}
			//this->m_Parameters[contour_id].bias[idx] = Components * log( 2*vnl_math::pi ) + log( det );

		} else {
			this->m_Parameters[contour_id].iCovariance[idx](0,0)=1.0/cov(0,0);
			//this->m_Parameters[contour_id].bias[idx] = Components * log( 2*vnl_math::pi ) + log( cov(0,0) );
		}
		this->m_Parameters[contour_id].bias[idx] = 1.0;
	}
}

template <typename TReferenceImageType, typename TCoordRepType>
typename MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::DeformationFieldPointer
MahalanobisLevelSets<TReferenceImageType,TCoordRepType>
::GetLevelSetsMap( MahalanobisLevelSets<TReferenceImageType,TCoordRepType>::DeformationFieldType* levelSetMap) {
	// Initialize interpolators
	this->m_Interp->SetInputImage( this->m_ReferenceImage );
	this->m_SparseToDenseResampler->CopyImageInformation( levelSetMap );


#ifndef NDEBUG
	double maxLS = -1.0;
#endif

	for( size_t cont = 0; cont < this->m_CurrentContourPosition.size(); cont++) {
		// Compute mesh of normals
		NormalFilterPointer normFilter = NormalFilterType::New();
		normFilter->SetInput( this->m_CurrentContourPosition[cont] );
		normFilter->Update();
		ContourDeformationPointer normals = normFilter->GetOutput();

		typename ContourDeformationType::PointsContainerConstIterator c_it = normals->GetPoints()->Begin();
		typename ContourDeformationType::PointsContainerConstIterator c_end = normals->GetPoints()->End();

		PointValueType sign[2] = { -1.0, 1.0 };
		PointValueType levelSet;
		PointType  ci;          //
		PointType  ci_prime;
		PixelType  fi;          // Feature on ci_prime
		typename ContourDeformationType::PixelType ni;
		typename ContourDeformationType::PointIdentifier idx;

		// for all node in mesh
		while (c_it!=c_end) {
			levelSet = 0;
			idx = c_it.Index();
			ci_prime = c_it.Value();
			fi = this->m_Interp->Evaluate( ci_prime );    // Feature in c'_{i}
			for( size_t i = 0; i<2; i++) {                 // Compute on both sides of the levelset
				PixelType dist = fi - this->m_Parameters[cont].mean[i];
				// compute mahalanobis distance in position
				levelSet+= sign[i] * ( this->m_Parameters[cont].bias[i] + dot_product(dist.GetVnlVector(), this->m_Parameters[cont].iCovariance[i].GetVnlMatrix() * dist.GetVnlVector() ) );
			}
			assert( !std::isnan(levelSet) );
			// project to normal, updating transform
			normals->GetPointData( idx, &ni );         // Normal ni in point c'_i
			ni*=levelSet;
			normals->SetPointData( idx, ni );
			++c_it;
#ifndef NDEBUG
			if ( fabs( levelSet ) > fabs( maxLS ) )
				maxLS = levelSet;
#endif
		}
		this->m_SparseToDenseResampler->SetInput( cont, normals );
	}

	// Interpolate sparse velocity field to targetDeformation
	this->m_SparseToDenseResampler->Update();
	return this->m_SparseToDenseResampler->GetOutput();

}

}

#endif /* MAHALANOBISLEVELSETS_HXX_ */
