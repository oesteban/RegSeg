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

#include <vnl/algo/vnl_sparse_lu.h>

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
}


template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetPhysicalDomainInformation( const DomainBase* image ) {
	size_t N = 1.0;
	for( size_t i=0; i<Dimension; i++) {
		N*=this->m_ControlPointsSize[i];
	}
	this->SetNumberOfParameters( N );

	if( N == 0 ){
		itkExceptionMacro( << "ControlPointsSize must be set and valid to set parameters this way.")
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
	this->InitializeCoefficientImages();
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::CopyGridInformation( const DomainBase* image ) {
	this->m_ControlPointsSize      = image->GetLargestPossibleRegion().GetSize();
	this->m_ControlPointsOrigin    = image->GetOrigin();
	this->m_ControlPointsSpacing   = image->GetSpacing();
	this->m_ControlPointsDirection = image->GetDirection();
	this->SetNumberOfParameters( image->GetLargestPossibleRegion().GetNumberOfPixels() );
	this->InitializeCoefficientImages();
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::InitializeCoefficientImages() {
	for( size_t dim = 0; dim < Dimension; dim++ ) {
		this->m_CoefficientImages[dim] = CoefficientsImageType::New();
		this->m_CoefficientImages[dim]->SetRegions(   this->m_ControlPointsSize );
		this->m_CoefficientImages[dim]->SetOrigin(    this->m_ControlPointsOrigin );
		this->m_CoefficientImages[dim]->SetSpacing(   this->m_ControlPointsSpacing );
		this->m_CoefficientImages[dim]->SetDirection( this->m_ControlPointsDirection );
		this->m_CoefficientImages[dim]->Allocate();
		this->m_CoefficientImages[dim]->FillBuffer( 0.0 );
	}

	if ( this->m_NumberOfParameters != this->m_CoefficientImages[0]->GetLargestPossibleRegion().GetNumberOfPixels() ) {
		itkExceptionMacro( << "inconsistent number of parameters.")
	}
}



template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetNumberOfSamples( size_t n ) {
	this->m_NumberOfSamples = n;
	this->m_OffGridPos.resize(n);
	for( size_t dim = 0; dim<Dimension; dim++ ) {
		this->m_OffGridValue[dim] = DimensionVector(n);
		this->m_OffGridValue[dim].fill( 0.0 );
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetNumberOfParameters( size_t K ) {
	this->m_NumberOfParameters = K;
	for( size_t dim = 0; dim<Dimension; dim++ ) {
		this->m_OnGridValue[dim] = DimensionVector(K);
		this->m_CoeffDerivative[dim] = DimensionVector(K);
		this->m_Coeff[dim] = DimensionVector(K);

		this->m_OnGridValue[dim].fill( 0.0 );
		this->m_CoeffDerivative[dim].fill( 0.0 );
		this->m_Coeff[dim].fill( 0.0 );

		for (size_t dim2 = 0; dim2<Dimension; dim2++ ) {
			this->m_Jacobian[dim][dim2] = DimensionVector(K);
			this->m_Jacobian[dim][dim2].fill( 0.0 );
		}
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputePhi() {
	this->m_Phi = WeightsMatrix(this->m_NumberOfSamples, this->m_NumberOfParameters);

	ScalarType wi;
	PointType ci, uk;
	size_t row, col;
	VectorType r;

	// Walk the grid region
	for( row = 0; row < this->m_Phi.rows(); row++ ) {
		ci = this->m_OffGridPos[row];
		for ( col=0; col < this->m_Phi.cols(); col++) {
			this->m_CoefficientImages[0]->TransformIndexToPhysicalPoint( this->m_CoefficientImages[0]->ComputeIndex( row ), uk );
			r = uk - ci;
			wi = this->Evaluate( r );

			if (wi > 0.0) {
				this->m_Phi.put(row, col, wi);
			}
		}
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeS() {
	PointType uk, uj;

	ScalarType wi;
	size_t row, col;
	VectorType r;

	this->m_S = WeightsMatrix(this->m_NumberOfParameters,this->m_NumberOfParameters);

	for( row = 0; row < this->m_S.rows(); row++ ) {
		this->m_CoefficientImages[0]->TransformIndexToPhysicalPoint( this->m_CoefficientImages[0]->ComputeIndex( row ), uk );

		for ( col=0; col < this->m_S.cols(); col++) {
			this->m_CoefficientImages[0]->TransformIndexToPhysicalPoint( this->m_CoefficientImages[0]->ComputeIndex( col ), uj );
			r = uk - uj;
			wi = this->Evaluate( r );

			if (wi > 0.0) {
				this->m_S.put(row, col, wi);
			}
		}
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeCoeffDerivatives() {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputePhi();
	}

    for ( size_t i = 0; i < Dimension; i++ ) {
    	this->m_CoeffDerivative[i].fill(0.0);
    	this->m_Phi.pre_mult(this->m_OffGridValue[i], this->m_CoeffDerivative[i] );
    }
}


template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::Interpolate() {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputePhi();
	}

	for( size_t i = 0; i < Dimension; i++ ) {
		this->m_OffGridValue[i].fill(0.0);
		// Interpolate
		this->m_Phi.mult( this->m_Coeff[i], this->m_OffGridValue[i] );
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeOnGridValues() {
	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeS();
	}

	for( size_t i = 0; i < Dimension; i++ ) {
		this->m_OnGridValue[i].fill(0.0);
		// Interpolate
		this->m_S.mult( this->m_Coeff[i], this->m_OnGridValue[i] );
	}
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeCoefficients() {
	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeS();
	}

	vnl_sparse_lu* system = new vnl_sparse_lu( this->m_S, vnl_sparse_lu::estimate_condition);

	for( size_t i = 0; i < Dimension; i++ ) {
		this->m_Coeff[i].fill( 0.0 );
		this->m_Coeff[i] = system->solve( this->m_OnGridValue[i] );

		ScalarType* cbuffer = this->m_CoefficientImages[i]->GetBufferPointer();
		for( size_t k = 0; k<this->m_NumberOfParameters; k++) {
			*( cbuffer + k ) = this->m_Coeff[i][k];
		}
	}
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
	for( size_t dim = 0; dim < Dimension; dim++) {
		ci[dim] = this->m_OffGridValue[dim][id];

		//if( std::isnan( ci[dim] ) || std::isinf( ci[dim] )) {
		//	ci[dim] = 0;
		//}
	}

	return ci;
}

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetOnGridValue( const size_t id ) {
	VectorType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_OnGridValue[dim][id];
		if( std::isnan( gi[dim] )) gi[dim] = 0;
	}
	return gi;
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

//template< class TScalarType, unsigned int NDimensions >
//inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
//SparseMatrixTransform<TScalarType,NDimensions>
//::GetOnGridWeight( const size_t id ) {
//	VectorType gi;
//	for( size_t dim = 0; dim < Dimension; dim++){
//		gi[dim] = this->m_CoeffDerivative[dim][id];
//		if( std::isnan( gi[dim] )) gi[dim] = 0;
//	}
//	return gi;
//}
//

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetCoefficient( const size_t id ) {
	VectorType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_Coeff[dim][id];
		if( std::isnan( gi[dim] )) gi[dim] = 0;
	}
	return gi;
}

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetCoeffDerivative( const size_t id ) {
	VectorType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_CoeffDerivative[dim][id];
		if( std::isnan( gi[dim] )) gi[dim] = 0;
	}
	return gi;
}

template< class TScalarType, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalarType,NDimensions>
::SetOffGridValue( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		if( this->m_OffGridValue[dim][id] != pi[dim] ) {
			this->m_OffGridValue[dim][id] = pi[dim];
			changed = true;
		}
	}
	return changed;
}

template< class TScalarType, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalarType,NDimensions>
::SetOnGridValue( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		if( this->m_OnGridValue[dim][id] != pi[dim]) {
			this->m_OnGridValue[dim][id] = pi[dim];
			changed=true;
		}
	}
	return changed;
}


template< class TScalarType, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalarType,NDimensions>
::SetCoefficient( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;

		// TODO check for a minimum size (it was done outside).
		if( this->m_Coeff[dim][id] != pi[dim] ) {
			this->m_Coeff[dim][id] = pi[dim];
			changed = true;
		}
	}
	return changed;
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
