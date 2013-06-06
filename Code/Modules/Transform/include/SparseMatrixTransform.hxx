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
::SparseMatrixTransform() {
	m_N = 0;
	m_k = 0;
}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputePhi() {
	this->m_k = this->m_GridPoints.size();
	this->m_N = this->m_ControlPoints.size();

	this->m_Phi = WeightsMatrix(this->m_k, this->m_N);
	ScalarType dist, wi;
	PointType ci, pi;
	size_t row, col;

	// Walk the grid region
	for( row = 0; row < this->m_k; row++ ) {
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
	this->m_InvertPhi = WeghtsMatrix( this->m_N, this->m_k );


}

template< class TScalarType, unsigned int NDimensions >
void
SparseMatrixTransform<TScalarType,NDimensions>
::ComputeGridPoints() {
	// TODO check m_Phi and initializations
	// TODO check if ControlPointsData changed


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
	// TODO check m_Phi and initializations
	// TODO check if ControlPointsData changed

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
::AddControlPoint(typename SparseMatrixTransform<TScalarType,NDimensions>::PointType pi ){
	this->m_ControlPoints.push_back( pi );
}

template< class TScalarType, unsigned int NDimensions >
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::AddGridPoint(typename SparseMatrixTransform<TScalarType,NDimensions>::PointType pi ){
	this->m_GridPoints.push_back( pi );
}
/*
template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::PointType&
SparseMatrixTransform<TScalarType,NDimensions>
::GetControlPoint( const size_t id ) {
	PointType ci;
	for( size_t dim = 0; dim < Dimension; dim++) {
		ci[dim] = this->m_ControlPoints[dim][id];
		if( std::isnan( ci[dim] )) ci[dim] = 0;
	}

	return ci;
}

template< class TScalarType, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalarType,NDimensions>::PointType&
SparseMatrixTransform<TScalarType,NDimensions>
::GetGridPoint( const size_t id ) {
	PointType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_GridPoints[dim][id];
		if( std::isnan( gi[dim] )) gi[dim] = 0;
	}

	return gi;
}
*/
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
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::SetControlPointData( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		this->m_ControlPointsData[dim][id] = pi[dim];
	}
}

template< class TScalarType, unsigned int NDimensions >
inline void
SparseMatrixTransform<TScalarType,NDimensions>
::SetGridPointData( const size_t id, typename SparseMatrixTransform<TScalarType,NDimensions>::VectorType pi ) {
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		this->m_GridPointsData[dim][id] = pi[dim];
	}
}

}


#endif /* SPARSEMATRIXTRANSFORM_HXX_ */
