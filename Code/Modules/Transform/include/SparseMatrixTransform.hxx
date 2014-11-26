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
SparseMatrixTransform<TScalar,NDimensions>
::SparseMatrixTransform():
Superclass(),
m_NumberOfParameters(0),
m_UseImageOutput(false) {
	this->m_ControlGridSize.Fill(10);
	this->m_ControlGridOrigin.Fill(0.0);
	this->m_ControlGridSpacing.Fill(0.0);
	this->m_ControlGridDirection.SetIdentity();
	this->m_ControlGridDirectionInverse.SetIdentity();

	this->m_Threader = itk::MultiThreader::New();
	this->m_NumberOfThreads = this->m_Threader->GetNumberOfThreads();

	this->m_ControlGridIndexToPhysicalPoint.SetIdentity();
	this->m_ControlGridPhysicalPointToIndex.SetIdentity();

	for( size_t i = 0; i<Dimension; i++ ) {
		this->m_PointValues[i] = DimensionVector();
	}
}

template< class TScalar, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalar,NDimensions>::ScalarType
SparseMatrixTransform<TScalar,NDimensions>
::EvaluateKernel( const VectorType r, const size_t dim ) {
	ScalarType wi=1.0;
	for (size_t i = 0; i<Dimension; i++) {
		wi*= this->m_KernelFunction->Evaluate( r[i] / this->m_ControlGridSpacing[i] );
	}
	return wi;
}

template< class TScalar, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalar,NDimensions>::ScalarType
SparseMatrixTransform<TScalar,NDimensions>
::EvaluateDerivative( const VectorType r, const size_t dim ) {
	ScalarType wi=1.0;
	for (size_t i = 0; i<Dimension; i++) {
		if( dim == i )
			wi*= this->m_DerivativeKernel->Evaluate( r[i] / this->m_ControlGridSpacing[i] );
		else
			wi*= this->m_KernelFunction->Evaluate( r[i] / this->m_ControlGridSpacing[i] );
	}
	return wi;
}

template< class TScalar, unsigned int NDimensions >
inline size_t
SparseMatrixTransform<TScalar,NDimensions>
::ComputeRegionOfPoint(const PointType& point, VectorType& cvector, IndexType& start, IndexType& end, OffsetTableType offsetTable ) {
	IndexType center;
	VectorType cstart,cend;
	size_t num = 1;
	offsetTable[0] = num;

    for ( size_t k = 0; k <Dimension; k++ ) {
    	cvector[k] = point[k] - this->m_ControlGridOrigin[k];
    }
    cvector = this->m_ControlGridPhysicalPointToIndex * cvector;

    for ( size_t k = 0; k <Dimension; k++ ) {
    	cstart[k] = cvector[k]-2.0;
    	cend[k] = cvector[k]+2.0;

    	if (cend[k] < 0.0 || cstart[k] > (this->m_ControlGridSize[k]-1) ) {
    		return 0;
    	}

    	center[k] = floor( cvector[k] + 0.5 );  // Add 0.5 to get an actual round() result
    	start[k] = ceil( cstart[k] );
    	end[k] = floor( cend[k] );

    	if ( start[k] < 0 ) start[k] = 0;
    	if ( end[k] > (this->m_ControlGridSize[k]-1) ) end[k] = this->m_ControlGridSize[k]-1;

    	if( end[k]<start[k] ) {
    		return 0;
    	}

    	num*= (  end[k]-start[k] + 1 );
    	offsetTable[k+1] = num;
    }

    return num;
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::SetControlGridInformation( const DomainBase* image ) {
	this->m_ControlGridSize      = image->GetLargestPossibleRegion().GetSize();
	this->m_ControlGridOrigin    = image->GetOrigin();
	this->m_ControlGridSpacing   = image->GetSpacing();
	this->m_ControlGridDirection = image->GetDirection();
	this->InitializeCoefficientsImages();
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::InitializeCoefficientsImages() {
	// Compute pixel conversion matrices
	MatrixType scale;


	for ( size_t i = 0; i < Dimension; i++ ) {
		if( this->m_ControlGridSize[i] == 0 ) {
			itkExceptionMacro(<< " control points grid size is not set.");
		}

		PointType step;
		if ( this->m_ControlGridSpacing[i] == 0.0 ) {
			step[i] = (this->m_DomainExtent[1][i] - this->m_DomainExtent[0][i])/ (1.0*this->m_ControlGridSize[i]);
			this->m_ControlGridSpacing[i] = fabs(step[i]);
			this->m_ControlGridOrigin[i] = this->m_DomainExtent[0][i] + 0.5 * step[i];
		}
		scale[i][i] = this->m_ControlGridSpacing[i];
	}

	if ( vnl_determinant( this->m_ControlGridDirection.GetVnlMatrix() ) == 0.0 ) {
	  itkExceptionMacro(<< "Bad direction, determinant is 0. Direction is " << this->m_ControlGridDirection);
	}

	typedef vnl_matrix< ScalarType > InternalComputationMatrix;
	InternalComputationMatrix dir( Dimension, Dimension );
	typedef vnl_matrix< typename DirectionType::ValueType > DirectionVNLMatrix;
	vnl_copy< DirectionVNLMatrix, InternalComputationMatrix >( this->m_ControlGridDirection.GetVnlMatrix(), dir );
	this->m_ControlGridIndexToPhysicalPoint = MatrixType(dir) * scale;
	this->m_ControlGridPhysicalPointToIndex = this->m_ControlGridIndexToPhysicalPoint.GetInverse();


	for( size_t dim = 0; dim < Dimension; dim++ ) {
		this->m_CoefficientsImages[dim] = CoefficientsImageType::New();
		this->m_CoefficientsImages[dim]->SetRegions(   this->m_ControlGridSize );
		this->m_CoefficientsImages[dim]->SetOrigin(    this->m_ControlGridOrigin );
		this->m_CoefficientsImages[dim]->SetSpacing(   this->m_ControlGridSpacing );
		this->m_CoefficientsImages[dim]->SetDirection( this->m_ControlGridDirection );
		this->m_CoefficientsImages[dim]->Allocate();
		this->m_CoefficientsImages[dim]->FillBuffer( 0.0 );
	}

	this->m_NumberOfParameters = this->m_CoefficientsImages[0]->GetLargestPossibleRegion().GetNumberOfPixels();

	PointType p;
	CoeffImageConstPointer ref =  itkDynamicCastInDebugMode< const CoefficientsImageType* >(this->m_CoefficientsImages[0].GetPointer() );
	for( size_t i = 0; i < this->m_NumberOfParameters; i++ ) {
		ref->TransformIndexToPhysicalPoint( ref->ComputeIndex( i ), p );
		this->m_ParamLocations.push_back( p );
	}

	VectorType zerov; zerov.Fill( 0.0 );
	this->m_Field = FieldType::New();
	this->m_Field->SetRegions(   this->m_ControlGridSize );
	this->m_Field->SetOrigin(    this->m_ControlGridOrigin );
	this->m_Field->SetSpacing(   this->m_ControlGridSpacing );
	this->m_Field->SetDirection( this->m_ControlGridDirection );
	this->m_Field->Allocate();
	this->m_Field->FillBuffer( zerov );

	for( size_t i = 0; i<Dimension; i++ ) {
		this->m_Derivatives[i] = CoefficientsImageType::New();
		this->m_Derivatives[i]->SetRegions(   this->m_ControlGridSize );
		this->m_Derivatives[i]->SetOrigin(    this->m_ControlGridOrigin );
		this->m_Derivatives[i]->SetSpacing(   this->m_ControlGridSpacing );
		this->m_Derivatives[i]->SetDirection( this->m_ControlGridDirection );
		this->m_Derivatives[i]->Allocate();
		this->m_Derivatives[i]->FillBuffer( 0.0 );
	}
	this->Modified();
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::ComputeMatrix( WeightsMatrixType type, size_t dim ) {
	struct SMTStruct str;
	str.Transform = this;
	str.type = type;
	str.dim = dim;
	size_t nCols = this->m_ParamLocations.size();

	switch( type ) {
	case Self::PHI:
		str.vrows = &this->m_PointLocations;
		if ( this->m_PointLocations.size() != this->m_NumberOfPoints ) {
			itkExceptionMacro(<< "OffGrid positions are not initialized");
		}
		this->m_Phi = WeightsMatrix( this->m_NumberOfPoints , nCols );
		str.matrix = &this->m_Phi;
		break;
	case Self::S:
		this->m_S = WeightsMatrix( nCols, nCols );

		str.vrows = &this->m_ParamLocations;
		str.matrix = &this->m_S;
		break;
	case Self::SPRIME:
		this->m_SPrime[dim] = WeightsMatrix( nCols, nCols );
		str.vrows = &this->m_ParamLocations;
		str.matrix = &this->m_SPrime[dim];
		break;
	default:
		itkExceptionMacro(<< "Matrix computation not implemented" );
		break;
	}

	this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
	this->GetMultiThreader()->SetSingleMethod( this->ComputeThreaderCallback, &str );
	this->GetMultiThreader()->SingleMethodExecute();
	//this->AfterThreadedComputeMatrix( &str );
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::AfterThreadedComputeMatrix( SMTStruct* str ) {}

template< class TScalar, unsigned int NDimensions >
ITK_THREAD_RETURN_TYPE
SparseMatrixTransform<TScalar,NDimensions>
::ComputeThreaderCallback(void *arg) {
	SMTStruct *str;

	itk::ThreadIdType total, threadId, threadCount;
	threadId = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
	threadCount = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;
	str = (SMTStruct *)( ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

	MatrixSectionType splitSection;
	splitSection.matrix = str->matrix;
	splitSection.vrows = str->vrows;
	splitSection.dim = str->dim;
	total = str->Transform->SplitMatrixSection( threadId, threadCount, splitSection );

	if( threadId < total ) {
		if (str->type == Self::SPRIME )
			str->Transform->ThreadedComputeMatrix( splitSection, &Self::EvaluateDerivative, threadId );
		else
			str->Transform->ThreadedComputeMatrix( splitSection, &Self::EvaluateKernel, threadId );
	}

	return ITK_THREAD_RETURN_VALUE;
}

template< class TScalar, unsigned int NDimensions >
itk::ThreadIdType
SparseMatrixTransform<TScalar,NDimensions>
::SplitMatrixSection( itk::ThreadIdType i, itk::ThreadIdType num, MatrixSectionType& section ) {
	size_t range = section.vrows->size();

	unsigned int valuesPerThread = itk::Math::Ceil< unsigned int >(range / (double)num);
	unsigned int maxThreadIdUsed = itk::Math::Ceil< unsigned int >(range / (double)valuesPerThread) - 1;

	section.section_id = i;
	section.first_row = i * valuesPerThread;

	if( i < maxThreadIdUsed ) {
		section.num_rows = valuesPerThread;
	}
	if( i == maxThreadIdUsed ) {
		section.num_rows = range - i*valuesPerThread;
	}

	return maxThreadIdUsed + 1;
}


template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::ThreadedComputeMatrix( MatrixSectionType& section, FunctionalCallback func, itk::ThreadIdType threadId ) {
	size_t last = section.first_row + section.num_rows;
	itk::SizeValueType nRows = section.num_rows;
	PointsList vrows = *(section.vrows);
	size_t dim = section.dim;

	ScalarType wi;
	PointType ci, uk;
	size_t row, col, number_of_pixels;
	VectorType r,cindex;
	IndexType start, end, current;
	OffsetTableType rOffsetTable;

	vcl_vector< int > cols;
	vcl_vector< ScalarType > vals;
	cols.resize(0);
	vals.resize(0);

	CoeffImagePointer ref = this->m_CoefficientsImages[0];


	// Walk the grid region
	for ( row = section.first_row; row < last; row++ ) {
		cols.clear();
		vals.clear();
		ci = vrows[row];
		number_of_pixels = this->ComputeRegionOfPoint( ci, cindex, start, end, rOffsetTable );

		for( size_t rOffset = 0; rOffset<number_of_pixels; rOffset++) {
			Helper::ComputeIndex( start, rOffset, rOffsetTable, current );
			TransformHelper::TransformIndexToPhysicalPoint( this->m_ControlGridIndexToPhysicalPoint, this->m_ControlGridOrigin, current, uk);
			r = ci - uk;
			wi = (this->*func)(r, dim);

			if ( fabs(wi) > 1.0e-5) {
				col = ref->ComputeOffset( current );
				cols.push_back( col );
				vals.push_back( wi );
			}
		}

		if ( cols.size() > 0 ) {
			section.matrix->set_row( row, cols, vals );
		}
	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::Interpolate( const DimensionParametersContainer& coeff ) {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputeMatrix( Self::PHI );
	}
	for( size_t i = 0; i<Dimension; i++ ) {
		this->m_PointValues[i].set_size( this->m_NumberOfPoints );
		this->m_Phi.mult( coeff[i], this->m_PointValues[i] );
	}

	if( this->m_UseImageOutput ) {
		bool setVector;
		ScalarType val;
		VectorType v; v.Fill(0.0);

		FieldPointer newfield = FieldType::New();
		newfield->SetRegions( this->m_DisplacementField->GetLargestPossibleRegion().GetSize() );
		newfield->SetOrigin( this->m_DisplacementField->GetOrigin() );
		newfield->SetSpacing( this->m_DisplacementField->GetSpacing() );
		newfield->SetDirection( this->m_DisplacementField->GetDirection() );
		newfield->Allocate();
		newfield->FillBuffer( v );

		VectorType* obuf = this->m_DisplacementField->GetBufferPointer();
		VectorType* ibuf = newfield->GetBufferPointer();

		for( size_t row = 0; row<this->m_NumberOfPoints; row++ ) {
			v.Fill( 0.0 );
			setVector = false;
			for( size_t i = 0; i<Dimension; i++) {
				val = this->m_PointValues[i][row];

				if( fabs(val) > 1.0e-5 ){
					v[i] = val;
					setVector = true;
				}
			}
			if (setVector) {
				*( obuf + row ) = v;
				*( ibuf + row ) = -v;
			}
		}
		this->SetInverseDisplacementField( newfield );

	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::InvertField() {
	// Check m_Phi_inverse and initializations
	if( this->m_Phi_inverse.rows() == 0 || this->m_Phi_inverse.cols() == 0 ) {
		this->InvertPhi();
	}

	DimensionParametersContainer coeff;
	// Compute coefficients
	for( size_t i = 0; i<Dimension; i++ ) {
		coeff[i] = DimensionVector( this->m_NumberOfPoints );
		this->m_Phi_inverse.mult( this->m_PointValues[i], coeff[i] );
	}

	// TODO set coeff here or inside Interpolate( coeff )

	// Interpolation with new coefficients
	this->UpdateField( coeff );
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::UpdateField( const DimensionParametersContainer& coeff ) {
	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeMatrix( Self::S );
	}

	DimensionParametersContainer fieldValues;

	// Interpolate
	for( size_t i = 0; i<Dimension; i++)
		this->m_S.mult( coeff[i], fieldValues[i] );

	VectorType v;
	ScalarType val;
	this->m_Field->FillBuffer( v );
	VectorType* fbuf = this->m_Field->GetBufferPointer();

	bool setVector;

	for ( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		v.Fill( 0.0 );
		setVector = false;
		for ( size_t col = 0; col < Dimension; col++ ) {
			val = fieldValues[col][row];
			if ( fabs(val) > 1.0e-5 ){
				v[col] = val;
				setVector = true;
			}
		}
		if (setVector)
			*( fbuf + row ) = v;
	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::ComputeCoefficients() {
	if ( this->m_NumberOfParameters == 0 ) {
		this->InitializeCoefficientsImages();
	}

	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeMatrix( Self::S );
	}

	SolverVector X[Dimension], Y[Dimension];

	DimensionParametersContainer fieldValues = this->VectorizeField( this->m_Field );
	DimensionParametersContainer coeffs;

	for ( size_t i = 0; i<Dimension; i++) {
		// TODO: check m_NumberOfPoints here
		// Y[i] = SolverVector( this->m_NumberOfPoints );
		Y[i] = SolverVector( fieldValues[i].size() );
		vnl_copy< DimensionVector, SolverVector >( fieldValues[i], Y[i] );
		X[i] = SolverVector( this->m_NumberOfParameters );
	}

	size_t nRows = this->m_S.rows();
	SolverMatrix S( nRows, this->m_S.cols() );
	SparseMatrixRowType row;
	typedef typename SolverMatrix::row SolverRow;
	SolverRow row_s;
	vcl_vector< int > cols;
	vcl_vector< double > vals;

	for( size_t i = 0; i < nRows; i++ ){
		row_s.clear();
		row = this->m_S.get_row( i );

		for( size_t j = 0; j< row.size(); j++ ) {
			cols.push_back( row[j].first );
			vals.push_back( static_cast< double >( row[j].second ) );

		}
		S.set_row( i, cols, vals );
	}

	if ( Dimension == 3 ) {
		SolverTypeTraits::Solve( S, Y[0], Y[1], Y[2], X[0], X[1], X[2] );
	} else if (Dimension == 2 ) {
		SolverTypeTraits::Solve( S, Y[0], Y[1], X[0], X[1] );
	} else {
		for( size_t col = 0; col < Dimension; col++ ) {
			SolverTypeTraits::Solve( S, Y[col], X[col] );
		}
	}

	for( size_t col = 0; col < Dimension; col++ ) {
		ScalarType* cbuffer = this->m_CoefficientsImages[col]->GetBufferPointer();
		for( size_t k = 0; k<this->m_NumberOfParameters; k++) {
			*( cbuffer + k ) = static_cast<ScalarType>(X[col][k]);
		}
	}

	this->Modified();
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::ComputeGradientField( ) {
	for( size_t i = 0; i<Dimension; i++ ) {
		if( this->m_SPrime[i].rows() == 0 || this->m_SPrime[i].cols() == 0 ) {
			this->ComputeMatrix( Self::SPRIME, i );
		}
	}

	DimensionParametersContainer coeff = this->VectorizeCoefficients();
	DimensionParametersContainer result[Dimension];
	ScalarType* fbuf[Dimension];


	for( size_t i = 0; i<Dimension; i++ ) {
		for( size_t j=0; j<Dimension; j++ ) {
			// Interpolate
			this->m_SPrime[j].mult( coeff[i], result[i][j] );
		}

		// Clear data buffer and get pointer
		this->m_Derivatives[i]->FillBuffer( 0.0 );
		fbuf[i] = this->m_Derivatives[i]->GetBufferPointer();
	}

	if (this->m_GradientField.IsNull()) {
		VectorType zero;
		zero.Fill(0.0);
		this->m_GradientField = FieldType::New();
		this->m_GradientField->SetRegions(this->m_Derivatives[0]->GetLargestPossibleRegion());
		this->m_GradientField->Allocate();
		this->m_GradientField->FillBuffer(zero);
	}
	VectorType* gbuf = this->m_GradientField->GetBufferPointer();


	VectorType v;
	VectorType g;
	ScalarType norm;

	for ( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
		g.Fill(0.0);
		for( size_t i = 0; i < Dimension; i++ ){
			v.Fill( 0.0 );

			for( size_t j = 0; j< Dimension; j++ ) {
				v[j] = result[i][j][row];
			}

			norm = v.GetSquaredNorm();
			if ( norm > 1.0e-7 )
				*( fbuf[i] + row ) = norm;
				g[i] = norm;
		}
		*( gbuf + row ) = g;
	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::InvertPhi() {
	// Check m_Phi
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputeMatrix( Self::PHI );
	}

	size_t nRows = this->m_Phi.rows();
	size_t nCols = this->m_Phi.cols();

	this->m_Phi_inverse = WeightsMatrix( nCols, nRows );

	ScalarType val;
	SolverMatrix A( nRows, nCols );
	SparseMatrixRowType row;
	typedef typename SolverMatrix::row SolverRow;
	SolverRow row_s;
	vcl_vector< int > cols;
	vcl_vector< double > vals;

	for( size_t i = 0; i < nRows; i++ ){
		row_s.clear();
		row = this->m_Phi.get_row( i );

		for( size_t j = 0; j< row.size(); j++ ) {
			cols.push_back( row[j].first );
			vals.push_back( static_cast< double >( row[j].second ) );
		}
		A.set_row( i, cols, vals );
	}

	SolverVector X, B;
	for( size_t col = 0; col<nRows; col++ ) {
		B = SolverVector( nCols );
		B.fill( 0.0 );
		B(col) = 1.0;
		X = SolverVector( nCols );

		SolverTypeTraits::Solve( A, B, X );

		for( size_t row = 0; row<nCols; row++ ) {
			val = X[row];
			if( fabs(val) > 1.0e-8 ) {
				this->m_Phi_inverse.put( row, col, val );
			}
		}
	}

	this->Modified();

}

template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::WeightsMatrix
SparseMatrixTransform<TScalar,NDimensions>
::GetPhi() {
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputeMatrix( Self::PHI );
	}

	return this->m_Phi;
}


//template< class TScalar, unsigned int NDimensions >
//void
//SparseMatrixTransform<TScalar,NDimensions>
//::SetParameters( const ParametersType & parameters ) {
//	// Save parameters. Needed for proper operation of TransformUpdateParameters.
//	//if (&parameters != &(this->m_Parameters)) {
//	//	this->m_Parameters = parameters;
//	//}
//    //
//	//if( this->m_PointLocations.size() == 0 ) {
//	//	itkExceptionMacro( << "Control Points should be set first." );
//	//}
//    //
//	//if( this->m_NumberOfPoints != this->m_PointLocations.size() ) {
//	//	if (! this->m_NumberOfPoints == 0 ){
//	//		itkExceptionMacro( << "N and number of control points should match." );
//	//	} else {
//	//		this->m_NumberOfPoints = this->m_PointLocations.size();
//	//	}
//	//}
//    //
//	//if ( this->m_NumberOfPoints != (parameters.Size() * Dimension) ) {
//	//	itkExceptionMacro( << "N and number of parameters should match." );
//	//}
//    //
//	//for ( size_t dim = 0; dim<Dimension; dim++) {
//	//	if ( this->m_PointValue[dim].size() == 0 ) {
//	//		this->m_PointValue[dim]( this->m_NumberOfPoints );
//	//	}
//	//	else if ( this->m_PointValue[dim].size()!=this->m_NumberOfPoints ) {
//	//		itkExceptionMacro( << "N and number of slots for parameter data should match." );
//	//	}
//	//}
//    //
//	//size_t didx = 0;
//    //
//	//while( didx < this->m_NumberOfPoints ) {
//	//	for ( size_t dim = 0; dim< Dimension; dim++ ) {
//	//		this->m_PointValue[dim][didx] = parameters[didx+dim];
//	//	}
//	//	didx++;
//	//}
//
//	// Modified is always called since we just have a pointer to the
//	// parameters and cannot know if the parameters have changed.
//	this->Modified();
//}



//template< class TScalar, unsigned int NDimensions >
//inline typename SparseMatrixTransform<TScalar,NDimensions>::JacobianType
//SparseMatrixTransform<TScalar,NDimensions>
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

template< class TScalar, unsigned int NDimensions >
inline typename SparseMatrixTransform<TScalar,NDimensions>::VectorType
SparseMatrixTransform<TScalar,NDimensions>
::GetCoefficient( const size_t id ) {
	VectorType gi;
	for( size_t dim = 0; dim < Dimension; dim++){
		gi[dim] = this->m_CoefficientsImages[dim]->GetPixel( this->m_CoefficientsImages[dim]->ComputeIndex( id ) );
	}
	return gi;
}

//template< class TScalar, unsigned int NDimensions >
//inline typename SparseMatrixTransform<TScalar,NDimensions>::VectorType
//SparseMatrixTransform<TScalar,NDimensions>
//::GetCoeffDerivative( const size_t id ) {
//	VectorType gi;
//	for( size_t dim = 0; dim < Dimension; dim++){
//		gi[dim] = this->m_CoeffDerivative[dim][id];
//		if( std::isnan( gi[dim] )) gi[dim] = 0;
//	}
//	return gi;
//}

template< class TScalar, unsigned int NDimensions >
inline bool
SparseMatrixTransform<TScalar,NDimensions>
::SetPointValue( const size_t id, typename SparseMatrixTransform<TScalar,NDimensions>::VectorType pi ) {
	bool changed = false;
	for( size_t dim = 0; dim < Dimension; dim++) {
		if( std::isnan( pi[dim] )) pi[dim] = 0;
		if( this->m_PointValues[dim][id] != pi[dim] ) {
			this->m_PointValues[dim][id] = pi[dim];
			changed = true;
		}
	}
	return changed;
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::SetCoefficientsImages( const CoefficientsImageArray & images ) {
	for( size_t dim = 0; dim < Dimension; dim++ ) {
		CoeffImageConstPointer c = itkDynamicCastInDebugMode< const CoefficientsImageType* >( images[dim].GetPointer() );
		this->SetCoefficientsImage( dim, c );
	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::SetCoefficientsVectorImage( const FieldType* images ) {
	ScalarType* buff[Dimension];

	if (this->m_NumberOfParameters == 0) {
		this->m_ControlGridDirection = images->GetDirection();
		this->m_ControlGridOrigin = images->GetOrigin();
		this->m_ControlGridSize = images->GetLargestPossibleRegion().GetSize();
		this->m_ControlGridSpacing = images->GetSpacing();
		this->InitializeCoefficientsImages();
	}

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

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::AddCoefficientsVectorImage( const FieldType* images ) {
	ScalarType* buff[Dimension];

	if (this->m_NumberOfParameters == 0) {
		this->SetCoefficientsVectorImage(images);
	} else {
		this->ComputeCoefficients();

		for( size_t dim = 0; dim < Dimension; dim++ ) {
			buff[dim] = this->m_CoefficientsImages[dim]->GetBufferPointer();
		}

		const VectorType* fbuf = images->GetBufferPointer();
		VectorType v;
		for( size_t row = 0; row < this->m_NumberOfParameters; row++ ) {
			v = *( fbuf + row );
			for( size_t dim = 0; dim < Dimension; dim++ ) {
				*( buff[dim] + row ) += v[dim];
			}
		}
	}
}

template< class TScalar, unsigned int NDimensions >
const typename SparseMatrixTransform<TScalar,NDimensions>::FieldType*
SparseMatrixTransform<TScalar,NDimensions>
::GetCoefficientsVectorImage(){
	if (this->m_CoefficientsField.IsNull()){
		this->m_CoefficientsField = FieldType::New();
		this->m_CoefficientsField->SetRegions(this->m_ControlGridSize);
		this->m_CoefficientsField->SetSpacing(this->m_ControlGridSpacing);
		this->m_CoefficientsField->SetDirection(this->m_ControlGridDirection);
		this->m_CoefficientsField->SetOrigin(this->m_ControlGridOrigin);
		this->m_CoefficientsField->Allocate();
	}

	const ScalarType* buff[Dimension];
	for( size_t dim = 0; dim < Dimension; dim++ ) {
		buff[dim] = this->m_CoefficientsImages[dim]->GetBufferPointer();
	}

	VectorType* fbuf = this->m_CoefficientsField->GetBufferPointer();
	VectorType v;
	for( size_t row = 0; row < this->m_NumberOfParameters; row++ ) {
		for( size_t dim = 0; dim < Dimension; dim++ ) {
			v[dim] = *( buff[dim] + row );
		}
		*(fbuf + row) = v;
	}
	return this->m_CoefficientsField;
}


template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
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

//template< class TScalar, unsigned int NDimensions >
//typename SparseMatrixTransform<TScalar,NDimensions>::WeightsMatrix
//SparseMatrixTransform<TScalar,NDimensions>
//::VectorizeCoefficients() {
//	WeightsMatrix m ( this->m_NumberOfParameters, Dimension );
//
//	std::vector< const ScalarType *> cbuffer;
//
//	for( size_t col = 0; col<Dimension; col++)
//		cbuffer.push_back( this->m_CoefficientsImages[col]->GetBufferPointer() );
//
//	ScalarType value;
//	for( size_t row = 0; row<this->m_NumberOfParameters; row++ ) {
//		for( size_t col = 0; col<Dimension; col++ ) {
//			value = *( cbuffer[col] + row );
//			if (value!=0) {
//				m.put( row, col, value );
//			}
//		}
//
//	}
//	return m;
//}


template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionParametersContainer
SparseMatrixTransform<TScalar,NDimensions>
::VectorizeCoefficients() const {
	DimensionParametersContainer m;
	const ScalarType* cbuffer;
	for( size_t i = 0; i<Dimension; i++) {
		m[i] = DimensionVector( this->m_NumberOfParameters );
		cbuffer = this->m_CoefficientsImages[i]->GetBufferPointer();

		for( size_t j = 0; j < this->m_NumberOfParameters; j++ ){
			m[i][j] = *( cbuffer + j );
		}
	}
	return m;
}

template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionParametersContainer
SparseMatrixTransform<TScalar,NDimensions>
::VectorizeDerivatives() const {
	DimensionParametersContainer m;
	const ScalarType* cbuffer;
	for( size_t i = 0; i<Dimension; i++) {
		m[i] = DimensionVector( this->m_NumberOfParameters );
		cbuffer = this->m_Derivatives[i]->GetBufferPointer();

		for( size_t j = 0; j < this->m_NumberOfParameters; j++ ){
			m[i][j] = *( cbuffer + j );
		}
	}
	return m;
}

template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionVector
SparseMatrixTransform<TScalar,NDimensions>
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
typename SparseMatrixTransform<TScalar,NDimensions>::WeightsMatrix
SparseMatrixTransform<TScalar,NDimensions>
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

template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionParametersContainer
SparseMatrixTransform<TScalar,NDimensions>
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


//template< class TScalar, unsigned int NDimensions >
//void
//SparseMatrixTransform<TScalar,NDimensions>
//::ComputeCoeffDerivatives() {
//	// Check m_Phi and initializations
//	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
//		this->ComputePhi();
//	}
//
//    for ( size_t i = 0; i < Dimension; i++ ) {
//    	this->m_CoeffDerivative[i].fill(0.0);
//    	this->m_Phi.pre_mult(this->m_PointValue[i], this->m_CoeffDerivative[i] );
//    }
//}

}


#endif /* SPARSEMATRIXTRANSFORM_HXX_ */
