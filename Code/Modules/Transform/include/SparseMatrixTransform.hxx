// --------------------------------------------------------------------------------------
// File:          SparseMatrixTransform.hxx
// Date:          Jun 6, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
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

#ifndef SPARSEMATRIXTRANSFORM_HXX_
#define SPARSEMATRIXTRANSFORM_HXX_

#include "SparseMatrixTransform.h"
#include <itkGaussianKernelFunction.h>
#include <itkBSplineKernelFunction.h>
#include <itkBSplineDerivativeKernelFunction.h>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageAlgorithm.h>

#include <vnl/algo/vnl_sparse_lu.h>
#include <vcl_vector.h>

namespace rstk {

template< class TScalarType, unsigned int NDimensions >
SparseMatrixTransform<TScalarType,NDimensions>
::SparseMatrixTransform(): Superclass(NDimensions) {
	this->m_NumberOfSamples = 0;
	this->m_NumberOfParameters = 0;

	this->m_ControlPointsSize.Fill(0);
	this->m_ControlPointsOrigin.Fill(0.0);
	this->m_ControlPointsSpacing.Fill(1.0);
	this->m_ControlPointsDirection.SetIdentity();
	this->m_ControlPointsDirectionInverse.SetIdentity();


	this->m_GridDataChanged = false;
	this->m_ControlDataChanged = false;
	this->m_UseImageOutput = false;
}

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::ScalarType
SparseMatrixTransform<TScalarType,NDimensions>
::Evaluate( const VectorType r ) const {
	ScalarType wi=1.0;
	for (size_t i = 0; i<Dimension; i++) {
		wi*= this->m_KernelFunction->Evaluate( fabs( r[i] ) / this->m_ControlPointsSpacing[i] );
	}
	return wi;
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetPhysicalDomainInformation( const DomainBase* image ) {
	for( size_t i=0; i<Dimension; i++) {
		if( this->m_ControlPointsSize[i] < 4 ){
			itkExceptionMacro( << "ControlPointsSize must be set and valid (>3) to set parameters this way.")
		}
	}

	ContinuousIndexType o_idx;
	o_idx.Fill( -0.5 );

	ContinuousIndexType e_idx;
	for ( size_t dim=0; dim< Dimension; dim++ ) {
		e_idx[dim] = image->GetLargestPossibleRegion().GetSize()[dim] - 0.5;
	}

	PointType first;
	PointType last;

	image->TransformContinuousIndexToPhysicalPoint( o_idx, first );
	image->TransformContinuousIndexToPhysicalPoint( e_idx, last );

	PointType orig,step;
	typename CoefficientsImageType::SpacingType spacing;
	for( size_t dim = 0; dim < Dimension; dim++ ) {
		step[dim] = (last[dim]-first[dim])/(1.0*this->m_ControlPointsSize[dim]);
		this->m_ControlPointsSpacing[dim] = fabs(step[dim]);
		this->m_ControlPointsOrigin[dim] = first[dim] + 0.5 * step[dim];
	}

	this->m_ControlPointsDirection = image->GetDirection();

	this->InitializeCoefficientsImages();
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetOutputReference( const DomainBase* image ) {
	this->m_UseImageOutput = true;
	this->m_OutputField = FieldType::New();
	this->m_OutputField->SetRegions( image->GetLargestPossibleRegion().GetSize() );
	this->m_OutputField->SetOrigin( image->GetOrigin() );
	this->m_OutputField->SetSpacing( image->GetSpacing() );
	this->m_OutputField->SetDirection( image->GetDirection() );
	this->m_OutputField->Allocate();
	VectorType zerov; zerov.Fill( 0.0 );
	this->m_OutputField->FillBuffer( zerov );
	size_t nparam = image->GetLargestPossibleRegion().GetNumberOfPixels();

	// Initialize off-grid positions
	PointType p;
	for( size_t i = 0; i < nparam; i++ ) {
		image->TransformIndexToPhysicalPoint( image->ComputeIndex( i ), p );
		this->m_OffGridPos.push_back( p );
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::CopyGridInformation( const DomainBase* image ) {
	this->m_ControlPointsSize      = image->GetLargestPossibleRegion().GetSize();
	this->m_ControlPointsOrigin    = image->GetOrigin();
	this->m_ControlPointsSpacing   = image->GetSpacing();
	this->m_ControlPointsDirection = image->GetDirection();
	this->InitializeCoefficientsImages();
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::InitializeCoefficientsImages() {
	for( size_t dim = 0; dim < Dimension; dim++ ) {
		this->m_CoefficientsImages[dim] = CoefficientsImageType::New();
		this->m_CoefficientsImages[dim]->SetRegions(   this->m_ControlPointsSize );
		this->m_CoefficientsImages[dim]->SetOrigin(    this->m_ControlPointsOrigin );
		this->m_CoefficientsImages[dim]->SetSpacing(   this->m_ControlPointsSpacing );
		this->m_CoefficientsImages[dim]->SetDirection( this->m_ControlPointsDirection );
		this->m_CoefficientsImages[dim]->Allocate();
		this->m_CoefficientsImages[dim]->FillBuffer( 0.0 );
	}

	this->m_NumberOfParameters = this->m_CoefficientsImages[0]->GetLargestPossibleRegion().GetNumberOfPixels();

	PointType p;
	CoeffImageConstPointer ref =  itkDynamicCastInDebugMode< const CoefficientsImageType* >(this->m_CoefficientsImages[0].GetPointer() );
	for( size_t i = 0; i < this->m_NumberOfParameters; i++ ) {
		ref->TransformIndexToPhysicalPoint( ref->ComputeIndex( i ), p );
		this->m_OnGridPos.push_back( p );
	}

	VectorType zerov; zerov.Fill( 0.0 );
	this->m_Field = FieldType::New();
	this->m_Field->SetRegions(   this->m_ControlPointsSize );
	this->m_Field->SetOrigin(    this->m_ControlPointsOrigin );
	this->m_Field->SetSpacing(   this->m_ControlPointsSpacing );
	this->m_Field->SetDirection( this->m_ControlPointsDirection );
	this->m_Field->Allocate();
	this->m_Field->FillBuffer( zerov );

}



template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetNumberOfSamples( size_t n ) {
	this->m_OffGridPos.resize(n);
}

template< class TScalarType, unsigned int NDimensions >
typename SparseMatrixTransform<TScalarType,NDimensions>::WeightsMatrix
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeMatrix( PointsList vrows, PointsList vcols ) {
	size_t nRows = vrows.size();
	size_t nCols = vcols.size();

	WeightsMatrix phi( nRows, nCols );

	ScalarType wi;
	PointType ci, uk;
	size_t row, col;
	VectorType r;

	vcl_vector< int > cols;
	vcl_vector< ScalarType > vals;

	cols.resize(0);
	vals.resize(0);

	// Walk the grid region
	for( row = 0; row < vrows.size(); row++ ) {
		cols.clear();
		vals.clear();

		ci = vrows[row];
		for ( col=0; col < vcols.size(); col++) {
			uk = vcols[col];
			r = uk - ci;
			wi = this->Evaluate( r );

			if (wi > 1.0e-5) {
				cols.push_back( col );
				vals.push_back( wi );
			}
		}

		if ( cols.size() > 0 ) {
			phi.set_row( row, cols, vals );
		}
	}
	return phi;
}


template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputePhi() {
	this->m_NumberOfSamples = this->m_OffGridPos.size();
	this->m_OffGridValueMatrix = WeightsMatrix ( this->m_OffGridPos.size(), Dimension );
	this->m_Phi = this->ComputeMatrix( this->m_OffGridPos, this->m_OnGridPos );
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeS() {
	this->m_S = this->ComputeMatrix( this->m_OnGridPos, this->m_OnGridPos );
}

//template< class TScalarType, unsigned int NDimensions >
//void
//SparseMatrixTransform<TScalarType,NDimensions>
//::ComputeCoeffDerivatives() {
//	// Check m_Phi and initializations
//	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
//		this->ComputePhi();
//	}
//
//    for ( size_t i = 0; i < Dimension; i++ ) {
//    	this->m_CoeffDerivative[i].fill(0.0);
//    	this->m_Phi.pre_mult(this->m_OffGridValue[i], this->m_CoeffDerivative[i] );
//    }
//}


template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::Interpolate() {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputePhi();
	}

	WeightsMatrix coeff = this->VectorizeCoefficients();
	this->m_OffGridValueMatrix = WeightsMatrix( this->m_NumberOfSamples, Dimension );
	this->m_Phi.mult( coeff, this->m_OffGridValueMatrix );

	SparseVectorType r;

	if( this->m_UseImageOutput ) {
		this->m_OutputField->FillBuffer( 0.0 );
		VectorType* obuf = this->m_OutputField->GetBufferPointer();
		VectorType v;
		bool setVector;

		for( size_t row = 0; row<this->m_NumberOfSamples; row++ ) {
			v.Fill( 0.0 );
			r = this->m_OffGridValueMatrix.get_row( row );
			setVector = false;

			for( size_t i = 0; i<r.size(); i++) {
				v[r[i].first] = r[i].second;
				setVector =  setVector || ( r[i].second != 0.0 );
			}
			if (setVector)
				*( obuf + row ) = v;
		}
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::UpdateField() {
	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeS();
	}

	WeightsMatrix coeff = this->VectorizeCoefficients();
	WeightsMatrix fieldValues;

	// Interpolate
	this->m_S.mult( coeff, fieldValues );

	VectorType v;
	v.Fill( 0.0 );
	this->m_Field->FillBuffer( v );
	VectorType* fbuf = this->m_Field->GetBufferPointer();

	SparseVectorType r;
	bool setVector;

	for ( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		setVector = false;
		r = fieldValues.get_row( row );
		for( size_t i = 0; i< r.size(); i++ ) {
			v[r[i].first] = r[i].second;
			setVector = setVector || ( r[i].second!= 0.0 );
		}

		if ( setVector )
			*( fbuf + row ) = v;
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeCoefficients() {
	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeS();
	}

	DimensionParametersContainer fieldValues = this->VectorizeField( this->m_Field );
	DimensionParametersContainer coeffs;

	if ( Dimension == 3 ) {
		SolverType::Solve( this->m_S, fieldValues[0], fieldValues[1], fieldValues[2], coeffs[0], coeffs[1], coeffs[2]  );
	} else if (Dimension == 2 ) {
		SolverType::Solve( this->m_S, fieldValues[0], fieldValues[1], coeffs[0], coeffs[1]  );
	} else {
		for( size_t col = 0; col < Dimension; col++ ) {
			SolverType::Solve( this->m_S, fieldValues[col], coeffs[col] );
		}
	}

	for( size_t col = 0; col < Dimension; col++ ) {
		ScalarType* cbuffer = this->m_CoefficientsImages[col]->GetBufferPointer();
		for( size_t k = 0; k<this->m_NumberOfParameters; k++) {
			*( cbuffer + k ) = coeffs[col][k];
		}
	}

	this->Modified();
}

//template< class TScalarType, unsigned int NDimensions >
//void
//SparseMatrixTransform<TScalarType,NDimensions>
//::ComputeJacobian( ) {
//	if( this->m_SPrime[0].rows() == 0 || this->m_SPrime[0].cols() == 0 ) {
//		this->ComputeSPrime();
//	}
//
//	for( size_t i = 0; i < Dimension; i++ ) {
//		for( size_t j = 0; j < Dimension; j++ ) {
//			this->m_SPrime[i].mult( this->m_Coeff[j], this->m_Jacobian[i][j] );
//		}
//	}
//}

template< class TScalarType, unsigned int NDimensions >
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::SetOffGridPos(size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::PointType pi ){
	this->m_OffGridPos[id] = pi;
}

template< class TScalarType, unsigned int NDimensions >
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::AddOffGridPos(typename SparseMatrixTransform<TScalarType,NDimensions>::PointType pi ){
	this->m_OffGridPos.push_back( pi );
}


template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetParameters( const ParametersType & parameters ) {
	// Save parameters. Needed for proper operation of TransformUpdateParameters.
	//if (&parameters != &(this->m_Parameters)) {
	//	this->m_Parameters = parameters;
	//}
    //
	//if( this->m_OffGridPos.size() == 0 ) {
	//	itkExceptionMacro( << "Control Points should be set first." );
	//}
    //
	//if( this->m_NumberOfSamples != this->m_OffGridPos.size() ) {
	//	if (! this->m_NumberOfSamples == 0 ){
	//		itkExceptionMacro( << "N and number of control points should match." );
	//	} else {
	//		this->m_NumberOfSamples = this->m_OffGridPos.size();
	//	}
	//}
    //
	//if ( this->m_NumberOfSamples != (parameters.Size() * Dimension) ) {
	//	itkExceptionMacro( << "N and number of parameters should match." );
	//}
    //
	//for ( size_t dim = 0; dim<Dimension; dim++) {
	//	if ( this->m_OffGridValue[dim].size() == 0 ) {
	//		this->m_OffGridValue[dim]( this->m_NumberOfSamples );
	//	}
	//	else if ( this->m_OffGridValue[dim].size()!=this->m_NumberOfSamples ) {
	//		itkExceptionMacro( << "N and number of slots for parameter data should match." );
	//	}
	//}
    //
	//size_t didx = 0;
    //
	//while( didx < this->m_NumberOfSamples ) {
	//	for ( size_t dim = 0; dim< Dimension; dim++ ) {
	//		this->m_OffGridValue[dim][didx] = parameters[didx+dim];
	//	}
	//	didx++;
	//}

	// Modified is always called since we just have a pointer to the
	// parameters and cannot know if the parameters have changed.
	this->Modified();
}


template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetOffGridValue( const size_t id ) const {
	VectorType ci;
	ci.Fill( 0.0 );

	if ( !this->m_OffGridValueMatrix.empty_row( id ) ) {
		for( size_t d = 0; d < Dimension; d++) {
			ci[d] = this->m_OffGridValueMatrix( id, d );
		}
	}

	return ci;
}

//template< class TScalarType, unsigned int NDimensions >
//inline typename SparseMatrixTransform<TScalarType,NDimensions>::JacobianType
//SparseMatrixTransform<TScalarType,NDimensions>
//::GetJacobian( const size_t id ) {
//	JacobianType gi;
//	gi.Fill( 0.0 );
//
//	for( size_t i = 0; i < Dimension; i++){
//		for( size_t j = 0; j < Dimension; j++){
//			gi(i,j) = this->m_Jacobian[j][i][id]; // Attention to transposition
//			if( std::isnan( gi(i,j) )) gi(i,j) = 0;
//		}
//	}
//	return gi;
//}

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetCoefficient( const size_t id ) {
	VectorType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_CoefficientsImages[dim]->GetPixel( this->m_CoefficientsImages[dim]->ComputeIndex( id ) );
	}
	return gi;
}

//template< class TScalarType, unsigned int NDimensions >
//inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
//SparseMatrixTransform<TScalarType,NDimensions>
//::GetCoeffDerivative( const size_t id ) {
//	VectorType gi;
//	for( size_t dim = 0; dim < Dimension; dim++){
//		gi[dim] = this->m_CoeffDerivative[dim][id];
//		if( std::isnan( gi[dim] )) gi[dim] = 0;
//	}
//	return gi;
//}

template< class TScalarType, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalarType,NDimensions>
::SetOffGridValue( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		if( this->m_OffGridValueMatrix.get(id, dim) != pi[dim] ) {
			this->m_OffGridValueMatrix.put(id, dim, pi[dim] );
			changed = true;
		}
	}
	return changed;
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetCoefficientsImages( const CoefficientsImageArray & images ) {
	for( size_t dim = 0; dim < Dimension; dim++ ) {
		CoeffImageConstPointer c = itkDynamicCastInDebugMode< const CoefficientsImageType* >( images[dim].GetPointer() );
		this->SetCoefficientsImage( dim, c );
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetCoefficientsVectorImage( const FieldType* images ) {
	ScalarType* buff[Dimension];

	for( size_t dim = 0; dim < Dimension; dim++ ) {
		buff[dim] = this->m_CoefficientsImages[dim]->GetBufferPointer();
	}

	const VectorType* fbuf = images->GetBufferPointer();
	VectorType v;
	for( size_t row = 0; row < this->m_NumberOfParameters; row++ ) {
		v = *( fbuf + row );
		for( size_t dim = 0; dim < Dimension; dim++ ) {
			*( buff[dim] + row ) = v[dim];
		}
	}
}


template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetCoefficientsImage( size_t dim, const CoefficientsImageType* c ) {
	if ( dim >= Dimension || dim < 0 ) {
		itkExceptionMacro(<< "trying to set coefficients for undefined dimension");
	}
	itk::ImageAlgorithm::Copy<CoefficientsImageType,CoefficientsImageType>(
			c, this->m_CoefficientsImages[dim],
			c->GetLargestPossibleRegion(),
			this->m_CoefficientsImages[dim]->GetLargestPossibleRegion()
	);
}

template< class TScalarType, unsigned int NDimensions >
typename SparseMatrixTransform<TScalarType,NDimensions>::WeightsMatrix
SparseMatrixTransform<TScalarType,NDimensions>
::VectorizeCoefficients() {
	WeightsMatrix m ( this->m_NumberOfParameters, Dimension );

	std::vector< const ScalarType *> cbuffer;

	for( size_t col = 0; col<Dimension; col++)
		cbuffer.push_back( this->m_CoefficientsImages[col]->GetBufferPointer() );

	ScalarType value;
	for( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		for( size_t col = 0; col<Dimension; col++ ) {
			value = *( cbuffer[col] + row );
			if (value!=0) {
				m.put( row, col, value );
			}
		}

	}
	return m;
}

template< class TScalarType, unsigned int NDimensions >
typename SparseMatrixTransform<TScalarType,NDimensions>::DimensionVector
SparseMatrixTransform<TScalarType,NDimensions>
::Vectorize( const CoefficientsImageType* image ) {
	DimensionVector v( image->GetLargestPossibleRegion().GetNumberOfPixels() );
	v.fill(0.0);

	const ScalarType *cbuffer = image->GetBufferPointer();

	for( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		v[row] = *( cbuffer + row );
	}
	return v;
}

template< class TScalarType, unsigned int NDimensions >
typename SparseMatrixTransform<TScalarType,NDimensions>::WeightsMatrix
SparseMatrixTransform<TScalarType,NDimensions>
::MatrixField( const FieldType* image ) {
	WeightsMatrix vectorized( image->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension );

	VectorType v;
	std::vector< unsigned int > cols;
	std::vector< ScalarType > vals;

	const VectorType *cbuffer = image->GetBufferPointer();
	for( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		v = *( cbuffer + row );
		cols.clear();
		vals.clear();

		for( size_t col = 0; col<Dimension; col++) {
			if( v[col]!= 0.0 ){
				cols.push_back( col );
				vals.push_back( v[col] );
			}
		}
		if ( cols.size() > 0 ) {
			vectorized.set_row( row, cols, vals );
		}
	}
	return vectorized;
}

template< class TScalarType, unsigned int NDimensions >
typename SparseMatrixTransform<TScalarType,NDimensions>::DimensionParametersContainer
SparseMatrixTransform<TScalarType,NDimensions>
::VectorizeField( const FieldType* image ) {
	DimensionParametersContainer vectorized;

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

//template< class TScalarType, unsigned int NDimensions >
//void
//SparseMatrixTransform<TScalarType,NDimensions>
//::ComputeSPrime( ) {
//	PointType s, x;
//	ScalarType wi;
//	size_t row, col;
//	VectorType r;
//
//	for ( size_t dim = 0; dim < Dimension; dim++ ) {
//		this->m_SPrime[dim] = WeightsMatrix(this->m_NumberOfParameters,this->m_NumberOfParameters);
//	}
//
//	for( row = 0; row < this->m_NumberOfParameters; row++ ) {
//		s = this->m_OnGridPos[row];
//
//		for ( col=0; col < this->m_NumberOfParameters; col++) {
//			x = this->m_OnGridPos[col];
//			r = s - x;
//			for ( size_t dim = 0; dim < Dimension; dim++ ) {
//				wi=1.0;
//				for (size_t i = 0; i<Dimension; i++) {
//					if( dim == i ) {
//						wi*= this->m_KernelDerivativeFunction->Evaluate( fabs( r[i] ) / m_Sigma[i] );
//					} else {
//						wi*= this->m_KernelFunction->Evaluate( fabs( r[i] ) / m_Sigma[i] );
//					}
//
//					if( fabs(wi) < vnl_math::eps ) {
//						wi = 0.0;
//						break;
//					}
//				}
//
//				if ( fabs(wi) > 0.0) {
//					//this->m_SPrime[dim].put(row, col, this->m_KernelNorm * wi);
//					this->m_SPrime[dim].put(row, col, wi);
//				}
//			}
//		}
//	}
//}

}

#endif /* SPARSEMATRIXTRANSFORM_HXX_ */
