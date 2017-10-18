// This file is part of RegSeg
//
// Copyright 2014-2017, Oscar Esteban <code@oscaresteban.es>
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#ifndef CACHEDMATRIXTRANSFORM_HXX_
#define CACHEDMATRIXTRANSFORM_HXX_

#include "CachedMatrixTransform.h"
#include <itkGaussianKernelFunction.h>
#include <itkBSplineKernelFunction.h>
#include <itkBSplineDerivativeKernelFunction.h>
#include <itkProgressReporter.h>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageAlgorithm.h>

#include <vnl/algo/vnl_sparse_lu.h>
#include <vnl/vnl_copy.h>
#include <vnl/vnl_matrix.h>
#include <vcl_vector.h>

namespace rstk {

template< class TScalar, unsigned int NDimensions >
CachedMatrixTransform<TScalar,NDimensions>
::CachedMatrixTransform():
Superclass(),
m_NumberOfPoints(0),
m_UseImageOutput(false),
m_InterpolationMode(UNKNOWN) {
	for( size_t i = 0; i<Dimension; i++ ) {
		this->m_PointValues[i] = DimensionVector();
	}
	this->m_ReferenceSize.Fill(0);
	this->m_ReferenceSpacing.Fill(0.0);
	this->m_ReferenceOrigin.Fill(0.0);
	this->m_ReferenceDirection.Fill(0.0);
	PointType zero; zero.Fill(0.0);
	this->m_DomainExtent.Fill(zero);
}

template< class TScalar, unsigned int NDimensions >
inline typename CachedMatrixTransform<TScalar,NDimensions>::VectorType
CachedMatrixTransform<TScalar,NDimensions>
::GetPointValue( const size_t id ) const {
	VectorType ci;
	ci.Fill( 0.0 );

	for( size_t d = 0; d < Dimension; d++) {
		ci[d] = this->m_PointValues[d][id];
	}

	return ci;
}


template< class TScalar, unsigned int NDimensions >
void
CachedMatrixTransform<TScalar,NDimensions>
::SetDomainExtent( const DomainBase* image ) {
	ContinuousIndexType o_idx;
	o_idx.Fill( -0.5 );

	ContinuousIndexType e_idx;
	for ( size_t dim=0; dim< Dimension; dim++ ) {
		e_idx[dim] = image->GetLargestPossibleRegion().GetSize()[dim] - 0.5;
	}

	image->TransformContinuousIndexToPhysicalPoint( o_idx, m_DomainExtent[0] );
	image->TransformContinuousIndexToPhysicalPoint( e_idx, m_DomainExtent[1] );

	for (size_t i = 0; i<Dimension; i++) {
		if( m_DomainExtent[1][i] < m_DomainExtent[0][i] ) {
			double tmp = m_DomainExtent[0][i];
			m_DomainExtent[0][i] = m_DomainExtent[1][i];
			m_DomainExtent[1][i] = tmp;
		}
	}
}


template< class TScalar, unsigned int NDimensions >
void
CachedMatrixTransform<TScalar,NDimensions>
::SetOutputReference( const DomainBase* image ) {
	if( this->m_InterpolationMode == POINTS_MODE ) {
		itkExceptionMacro(<< "trying to change iterpolation mode from scattered to regular.");
	} else {
		this->m_InterpolationMode = GRID_MODE;
	}

	this->m_UseImageOutput = true;
	VectorType zerov; zerov.Fill( 0.0 );

	this->m_ReferenceSize = image->GetLargestPossibleRegion().GetSize();
	this->m_ReferenceSpacing = image->GetSpacing();
	this->m_ReferenceOrigin = image->GetOrigin();
	this->m_ReferenceDirection = image->GetDirection();
	SetFieldParametersFromImage(image);

	this->m_NumberOfPoints = image->GetLargestPossibleRegion().GetNumberOfPixels();
	// Initialize off-grid positions
	PointType p;
	for( size_t i = 0; i < this->m_NumberOfPoints; i++ ) {
		image->TransformIndexToPhysicalPoint( image->ComputeIndex( i ), p );
		this->m_PointLocations.push_back( p );
	}
}

template< class TScalar, unsigned int NDimensions >
void CachedMatrixTransform<TScalar,NDimensions>
::SetFieldParametersFromImage(const DomainBase* image) {
	SizeType s = image->GetLargestPossibleRegion().GetSize();
	PointType o = image->GetOrigin();
	SpacingType sp = image->GetSpacing();
	DirectionType d = image->GetDirection();

	ParametersType param;
	param.SetSize(Dimension * (Dimension + 3));

	for(size_t i = 0; i<Dimension; i++) {
		param[i] = s[i];
		param[i + Dimension] = o[i];
		param[i + Dimension * 2] = sp[i];
		param[i + Dimension * 3] = d[0][i];
		param[i + Dimension * 4] = d[1][i];
		param[i + Dimension * 5] = d[2][i];
	}
	this->SetFieldFixedParameters(param);
}

template< class TScalar, unsigned int NDimensions >
void CachedMatrixTransform<TScalar,NDimensions>
::SetCoefficientsParametersFromImage(const DomainBase* image) {
	SizeType s = image->GetLargestPossibleRegion().GetSize();
	PointType o = image->GetOrigin();
	SpacingType sp = image->GetSpacing();
	DirectionType d = image->GetDirection();

	ParametersType param;
	param.SetSize(Dimension * (Dimension + 3));

	for(size_t i = 0; i<Dimension; i++) {
		param[i] = s[i];
		param[i + Dimension] = o[i];
		param[i + Dimension * 2] = sp[i];
		param[i + Dimension * 3] = d[0][i];
		param[i + Dimension * 4] = d[1][i];
		param[i + Dimension * 5] = d[2][i];
	}
	this->SetCoefficientsFixedParameters(param);
}

template< class TScalar, unsigned int NDimensions >
void CachedMatrixTransform<TScalar,NDimensions>
::SetOutputPoints( const PointsList points ){
	if( this->m_InterpolationMode ==  GRID_MODE) {
		itkExceptionMacro(<< "trying to change iterpolation mode from regular to scattered.");
	} else {
		this->m_InterpolationMode = POINTS_MODE;
	}

	this->m_PointLocations = points;
	this->m_NumberOfPoints = points.size();
	this->Modified();
}

template< class TScalar, unsigned int NDimensions >
void CachedMatrixTransform<TScalar,NDimensions>
::SetOutputPoints( const PointsList points, const PointIdContainer valid ){
	this->SetOutputPoints(points);
	this->m_ValidLocations = valid;
}


template< class TScalar, unsigned int NDimensions >
typename CachedMatrixTransform<TScalar,NDimensions>::DimensionVector
CachedMatrixTransform<TScalar,NDimensions>
::Vectorize( const CoefficientsImageType* image ) {
	DimensionVector v( image->GetLargestPossibleRegion().GetNumberOfPixels() );
	v.fill(0.0);

	const ScalarType *cbuffer = image->GetBufferPointer();

	for( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		v[row] = *( cbuffer + row );
	}
	return v;
}

template< class TScalar, unsigned int NDimensions >
typename CachedMatrixTransform<TScalar,NDimensions>::DimensionParameters
CachedMatrixTransform<TScalar,NDimensions>
::VectorizeField( const FieldType* image ) {
	DimensionParameters vectorized;

	for( size_t col = 0; col<Dimension; col++) {
		vectorized[col] = DimensionVector( image->GetLargestPossibleRegion().GetNumberOfPixels() );
		vectorized[col].fill(0.0);
	}

	VectorType v;
	const VectorType *cbuffer = image->GetBufferPointer();

	for( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		v = *( cbuffer + row );
		for( size_t col = 0; col<Dimension; col++) {
			vectorized[col][row] = v[col];
		}
	}
	return vectorized;
}

}


#endif /* CACHEDMATRIXTRANSFORM_HXX_ */
