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

#ifndef SPARSETODENSEFIELDRESAMPLEFILTER_HXX_
#define SPARSETODENSEFIELDRESAMPLEFILTER_HXX_

#include "SparseToDenseFieldResampleFilter.h"
#include "IDWMultivariateInterpolator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace rstk {

template<class TInputMesh, class TDenseFieldType>
SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::SparseToDenseFieldResampleFilter() {
	m_FieldSpacing.Fill(1.0);
	m_FieldOrigin.Fill(0.0);
	m_FieldDirection.SetIdentity();
	m_FieldSize.Fill(0);
	m_FieldStartIndex.Fill(0);
	m_IsPhiInitialized = false;
	m_N = 0;
	m_k = 0;
};

template<class TInputMesh, class TDenseFieldType>
inline void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::
AddControlPoint( const typename SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::PointType p ) {
//	FieldContinuousIndexType idx;
//
//	if( ! this->m_Field->TransformPhysicalPointToContinuousIndex( p, idx ) ) {
//		itkExceptionMacro( << "ControlPoint p=[" << p << "] is outside deformation field extents." );
//	}
	this->m_ControlPoints.push_back( p );
	this->m_N+=1;
}

template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::
AddControlPoints( const typename SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::ContourType* prior) {
	for( size_t i = 0; i<prior->GetNumberOfPoints();i++ ) {
		this->AddControlPoint( prior->GetPoint(i) );
	}
}

template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::PrintSelf(
		std::ostream &os, Indent indent) const {
	//Superclass::PrintSelf(os, indent);

	os << indent << "Size: " << m_FieldSize << std::endl;
	os << indent << "FieldStartIndex: " << m_FieldStartIndex << std::endl;
	os << indent << "FieldSpacing: " << m_FieldSpacing << std::endl;
	os << indent << "FieldOrigin: " << m_FieldOrigin << std::endl;
	os << indent << "FieldDirection: " << m_FieldDirection << std::endl;
	//os << indent << "Interpolator: " << m_Interpolator.GetPointer()
	//		<< std::endl;
}
;

/**
 * Set the output image spacing.
 */
template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::SetFieldSpacing(
		const double *spacing) {
	FieldSpacingType s(spacing);

	this->SetFieldSpacing(s);
}

/**
 * Set the output image origin.
 */
template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::SetFieldOrigin(
		const double *origin) {
	FieldPointType p(origin);
	this->SetFieldOrigin(p);
}

/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::GenerateData() {
	if ( m_N == 0 ) {
		itkExceptionMacro(<< "Shape priors not set");
	}

	FieldPointer field = this->GetOutput();

	m_k = field->GetLargestPossibleRegion().GetNumberOfPixels();

	for (size_t i = 0; i < Dimension; i++ ) m_LevelSetVector[i] = GradientsVector( m_N );
	size_t nConts = this->GetNumberOfIndexedInputs();

	if (!m_IsPhiInitialized) {
		m_Phi = WeightsMatrix(m_k, m_N);
		double dist, wi;
		FieldPointType ci, pi;
		size_t row, col;

		// Walk the output region
		for( row = 0; row < this->m_k; row++ ) {
			field->TransformIndexToPhysicalPoint(field->ComputeIndex(row), pi);

			for ( col=0; col < this->m_N; col++) {
				// TODO Extract RBF from here (use a functor?)
				// Compute weight i: wi(x) (inverse distance)
				dist = (pi - this->m_ControlPoints[col]).GetNorm();
				wi = (dist < vnl_math::eps) ? 1.0 : (1.0 / vcl_pow(dist, 2));
				if (wi > 1e-3) {
					m_Phi.put(row, col, wi);
				}
			}
		}

		m_IsPhiInitialized = true;
	}

	size_t contsOffset = 0;
	VectorType ni;
	for( size_t cont = 0; cont< nConts; cont++ ) {
		size_t nPoints = this->GetInput(cont)->GetNumberOfPoints();
		for (size_t el = 0; el < nPoints; el++) {
			this->GetInput(cont)->GetPointData( el, &ni );

			for ( size_t i = 0; i < Dimension; i++ ) {
				if( std::isnan( ni[i] )) ni[i] = 0;
				m_LevelSetVector[i].put( contsOffset + el, ni[i] );
			}
		}
		contsOffset += nPoints;
	}

	double norms = 0.0;
	for (size_t i = 0; i < Dimension; i++ ) {
		norms+= m_LevelSetVector[i].one_norm();
	}


	if( norms==0 ) return;

    GradientsVector outvector[3];
    norms=0;
    for ( size_t i = 0; i < Dimension; i++ ) {
    	m_Phi.mult(m_LevelSetVector[i], outvector[i] );
    	norms+= outvector[i].one_norm();
    }

    if( norms==0 ) return;


    VectorType* buffer = field->GetBufferPointer();

    for( size_t k = 0; k<m_k; k++ ) {
    	ni = itk::NumericTraits<VectorType>::Zero;
    	for ( size_t i = 0; i < Dimension; i++ ) {
    		ni[i] = outvector[i].get(k);
    	}

    	if ( ni.GetNorm()!=0 ){
    		*( buffer + k ) = ni;
    	}
    }

}

template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>
::CopyImageInformation( const FieldType* image) {
	m_FieldDirection = image->GetDirection();
	m_FieldOrigin = image->GetOrigin();
	m_FieldSpacing = image->GetSpacing();
	m_FieldStartIndex = image->GetRequestedRegion().GetIndex();
	m_FieldSize = image->GetRequestedRegion().GetSize();
}

/**
 * Inform pipeline of required output region
 */
template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>
::UpdateOutputInformation() {
	this->Modified();
	this->Superclass::UpdateOutputInformation();
}

/**
 * Inform pipeline of required output region
 */
template<class TInputMesh, class TDenseFieldType>
void SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>
::GenerateOutputInformation() {
	// call the superclass' implementation of this method
	this->Superclass::GenerateOutputInformation();

	// get pointers to the input and output
	m_Field=this->GetOutput();

	if (!m_Field) {
		return;
	}

	// Set the size of the output region
	typename FieldType::RegionType outputLargestPossibleRegion;
	outputLargestPossibleRegion.SetSize(m_FieldSize);
	outputLargestPossibleRegion.SetIndex(m_FieldStartIndex);
	m_Field->SetLargestPossibleRegion(outputLargestPossibleRegion);
	m_Field->SetBufferedRegion(outputLargestPossibleRegion);
	m_Field->Allocate();
	m_Field->FillBuffer( itk::NumericTraits<VectorType>::Zero );

	// Set spacing and origin
	m_Field->SetSpacing(m_FieldSpacing);
	m_Field->SetOrigin(m_FieldOrigin);
	m_Field->SetDirection(m_FieldDirection);

	return;
}

/**
 * Verify if any of the components has been modified.
 */
template<class TInputMesh, class TDenseFieldType>
unsigned long SparseToDenseFieldResampleFilter<TInputMesh, TDenseFieldType>::GetMTime(
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
