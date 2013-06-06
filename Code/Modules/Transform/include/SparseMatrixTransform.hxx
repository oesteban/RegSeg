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

namespace rstk {

template< class TScalarType, unsigned int NDimensions >
SparseMatrixTransform<TScalarType,NDimensions>
::SparseMatrixTransform(): Superclass(NDimensions) {
	this->m_N = 0;
	this->m_K = 0;
	this->m_GridDataChanged = false;
	this->m_ControlDataChanged = false;
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetN( size_t N ) {
	this->m_N = N;
	this->m_ControlPoints.resize(N);
	for( size_t dim = 0; dim<Dimension; dim++ )
		this->m_ControlPointsData[dim] = DimensionVector(N);

}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::SetK( size_t K ) {
	this->m_K = K;
	this->m_GridPoints.resize(K);
	for( size_t dim = 0; dim<Dimension; dim++ )
		this->m_GridPointsData[dim] = DimensionVector(K);

	for ( size_t k = 0; k<K; k++ ) {
		this->SetGridPointData( k, itk::NumericTraits<VectorType>::Zero );
	}

}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputePhi() {
	this->m_Phi = WeightsMatrix(this->m_K, this->m_N);
	this->m_InvertPhi = WeightsMatrix( this->m_N, this->m_K );

	ScalarType dist, wi;
	PointType ci, pi;
	size_t row, col;

	// Walk the grid region
	for( row = 0; row < this->m_K; row++ ) {
		pi = this->m_GridPoints[row];

		for ( col=0; col < this->m_N; col++) {
			ci = this->m_ControlPoints[col];
			wi = this->ComputeWeight( pi, ci );
			if (wi > 1e-3) {
				this->m_Phi.put(row, col, wi);
			}
		}
	}

	// TODO Invert phi matrix
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeGridPoints() {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputePhi();
	}

	// Check if ControlPointsData changed
	double norms = 0.0;
	for (size_t i = 0; i < Dimension; i++ ) {
		norms+= this->m_ControlPointsData[i].one_norm();
	}

	if( norms==0 ) return; // Nothing to do

    for ( size_t i = 0; i < Dimension; i++ ) {
    	this->m_Phi.mult(this->m_ControlPointsData[i], this->m_GridPointsData[i] );
    }
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeControlPoints() {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputePhi();
	}

	// Check if GridPointsData changed
	double norms = 0.0;
	for (size_t i = 0; i < Dimension; i++ ) {
		norms+= this->m_GridPointsData[i].one_norm();
	}

	if( norms==0 ) return; // Nothing to do

    for ( size_t i = 0; i < Dimension; i++ ) {
    	this->m_InvertPhi.mult(this->m_GridPointsData[i], this->m_ControlPointsData[i] );
    }

}


template< class TScalarType, unsigned int NDimensions >
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::SetControlPoint(size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::PointType pi ){
	this->m_ControlPoints[id] = pi;
}

template< class TScalarType, unsigned int NDimensions >
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::SetGridPoint(size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::PointType pi ){
	this->m_GridPoints[id] = pi;
}

template<class TScalarType, unsigned int NDimensions>
void SparseMatrixTransform<TScalarType, NDimensions>::SetParameters(
		const ParametersType & parameters) {
	// Save parameters. Needed for proper operation of TransformUpdateParameters.
	if (&parameters != &(this->m_Parameters)) {
		this->m_Parameters = parameters;
	}

	if( this->m_ControlPoints.size() == 0 ) {
		itkExceptionMacro( << "Control Points should be set first." );
	}

	if( this->m_N != this->m_ControlPoints.size() ) {
		if (! this->m_N == 0 ){
			itkExceptionMacro( << "N and number of control points should match." );
		} else {
			this->m_N = this->m_ControlPoints.size();
		}
	}

	if ( this->m_N != (parameters.Size() * Dimension) ) {
		itkExceptionMacro( << "N and number of parameters should match." );
	}

	for ( size_t dim = 0; dim<Dimension; dim++) {
		if ( this->m_ControlPointsData[dim].size() == 0 ) {
			this->m_ControlPointsData[dim]( this->m_N );
		}
		else if ( this->m_ControlPointsData[dim].size()!=this->m_N ) {
			itkExceptionMacro( << "N and number of slots for parameter data should match." );
		}
	}

	size_t didx = 0;

	while( didx < this->m_N ) {
		for ( size_t dim = 0; dim< Dimension; dim++ ) {
			this->m_ControlPointsData[dim][didx] = parameters[didx+dim];
		}
		didx++;
	}

	// Modified is always called since we just have a pointer to the
	// parameters and cannot know if the parameters have changed.
	this->Modified();
}


template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetControlPointData( const size_t id ) {
	VectorType ci;
	for( size_t dim = 0; dim < Dimension; dim++) {
		ci[dim] = this->m_ControlPointsData[dim][id];
		if( std::isnan( ci[dim] )) ci[dim] = 0;
	}

	return ci;
}

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType
SparseMatrixTransform<TScalarType,NDimensions>
::GetGridPointData( const size_t id ) {
	VectorType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_GridPointsData[dim][id];
		if( std::isnan( gi[dim] )) gi[dim] = 0;
	}

	return gi;
}

template< class TScalarType, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalarType,NDimensions>
::SetControlPointData( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		if( this->m_ControlPointsData[dim][id] != pi[dim] ) {
			this->m_ControlPointsData[dim][id] = pi[dim];
			changed = true;
		}
	}
	return changed;
}

template< class TScalarType, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalarType,NDimensions>
::SetGridPointData( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		if( this->m_GridPointsData[dim][id] != pi[dim]) {
			this->m_GridPointsData[dim][id] = pi[dim];
			changed=true;
		}
	}
	return changed;
}

}


#endif /* SPARSEMATRIXTRANSFORM_HXX_ */
