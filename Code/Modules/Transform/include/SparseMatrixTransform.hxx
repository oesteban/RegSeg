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
// This file is part of RegSeg
//
// RegSeg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RegSeg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RegSeg.  If not, see <http://www.gnu.org/licenses/>.
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
Superclass() {
	this->m_ControlGridSize.Fill(10);
	this->m_ControlGridOrigin.Fill(0.0);
	this->m_ControlGridSpacing.Fill(0.0);
	this->m_ControlGridDirection.SetIdentity();
	this->m_ControlGridDirectionInverse.SetIdentity();
	this->m_MaximumDisplacement.Fill(0.0);

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
	this->InitializeCoefficientsImages();
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::InitializeCoefficientsImages() {
	// Compute pixel conversion matrices
	MatrixType scale;
	ParametersType param;

	for ( size_t i = 0; i < Dimension; i++ ) {
		if( this->m_ControlGridSize[i] == 0 ) {
			itkExceptionMacro(<< " control points grid size is not set.");
		}

		double extent = fabs(this->m_DomainExtent[1][i] - this->m_DomainExtent[0][i]);

		if ( extent <= 1e-3) {
			itkExceptionMacro(<< " extent is not defined.");
		}

		if ( this->m_ControlGridSpacing[i] == 0.0 ) {
			this->m_ControlGridSpacing[i] = extent / (1.0*this->m_ControlGridSize[i]);
			this->m_ControlGridOrigin[i] = this->m_DomainExtent[0][i] + 0.5 * this->m_ControlGridSpacing[i];
		}
		scale[i][i] = this->m_ControlGridSpacing[i];
	}

	if ( vnl_determinant( this->m_ControlGridDirection.GetVnlMatrix() ) == 0.0 ) {
	  itkExceptionMacro(<< "Bad direction, determinant is 0. Direction is " << this->m_ControlGridDirection);
	}

	param.SetSize(Dimension * (Dimension + 3));
	for(size_t i = 0; i<Dimension; i++) {
		param[i] = this->m_ControlGridSize[i];
		param[i + Dimension] = this->m_ControlGridOrigin[i];
		param[i + Dimension * 2] = this->m_ControlGridSpacing[i];
		param[i + Dimension * 3] = this->m_ControlGridDirection[0][i];
		param[i + Dimension * 4] = this->m_ControlGridDirection[1][i];
		param[i + Dimension * 5] = this->m_ControlGridDirection[2][i];
	}
	this->SetCoefficientsFixedParameters(param);

	typedef vnl_matrix< ScalarType > InternalComputationMatrix;
	InternalComputationMatrix dir( Dimension, Dimension );
	typedef vnl_matrix< typename DirectionType::ValueType > DirectionVNLMatrix;
	vnl_copy< DirectionVNLMatrix, InternalComputationMatrix >( this->m_ControlGridDirection.GetVnlMatrix(), dir );
	this->m_ControlGridIndexToPhysicalPoint = MatrixType(dir) * scale;
	this->m_ControlGridPhysicalPointToIndex = this->m_ControlGridIndexToPhysicalPoint.GetInverse();

	for( size_t i = 0; i<Dimension; i++ ) {
		this->m_Derivatives[i] = CoefficientsImageType::New();
		this->m_Derivatives[i]->SetRegions(   this->m_ControlGridSize );
		this->m_Derivatives[i]->SetOrigin(    this->m_ControlGridOrigin );
		this->m_Derivatives[i]->SetSpacing(   this->m_ControlGridSpacing );
		this->m_Derivatives[i]->SetDirection( this->m_ControlGridDirection );
		this->m_Derivatives[i]->Allocate();
		this->m_Derivatives[i]->FillBuffer( 0.0 );
	}

	for (size_t i = 0; i<Dimension; i++) {
		this->m_MaximumDisplacement[i] = 0.40 * this->m_ControlGridSpacing[i];
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
		str.vcols = &this->m_ParamLocations;
		if ( this->m_PointLocations.size() != this->m_NumberOfPoints ) {
			itkExceptionMacro(<< "OffGrid positions are not initialized");
		}
		this->m_Phi = WeightsMatrix( this->m_NumberOfPoints , nCols );
		str.matrix = &this->m_Phi;
		break;

	case Self::PHI_FIELD:
		str.vrows = &this->m_FieldLocations;
		str.vcols = &this->m_ParamLocations;
		this->m_FieldPhi = WeightsMatrix( this->m_FieldLocations.size() , nCols );
		str.matrix = &this->m_FieldPhi;
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
	this->AfterComputeMatrix(type);
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::AfterComputeMatrix(WeightsMatrixType type) {
	size_t numvalid = this->m_ValidLocations.size();
	this->m_Phi.normalize_rows();
	if(numvalid > 0 &&  numvalid <= this->m_NumberOfPoints ) {
		if (type == Self::PHI)  {
			this->m_Phi_valid = WeightsMatrix( numvalid , this->m_NumberOfDimParameters );

			size_t rid = 0;
			typename PointIdContainer::iterator it = this->m_ValidLocations.begin();
			typename PointIdContainer::iterator end = this->m_ValidLocations.end();
			typename WeightsMatrix::row row;
			vcl_vector< int > cols;
			vcl_vector< ScalarType > vals;

			while( it != end ) {
				row = this->m_Phi.get_row(*it);
				cols.resize(0);
				vals.resize(0);

				for(size_t i = 0; i < row.size(); i++) {
					cols.push_back(row[i].first);
					vals.push_back(row[i].second);
				}

				this->m_Phi_valid.set_row(rid, cols, vals);
				rid++;
				++it;
			}
		}
	}


}

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

	FieldPointer ref = this->m_CoefficientsField;
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
::InterpolatePoints() {
	const DimensionParameters coeff = this->VectorizeCoefficients();
	// Check m_Phi and initializations
	if( this->m_Phi.rows() == 0 || this->m_Phi.cols() == 0 ) {
		this->ComputeMatrix( Self::PHI );
	}
	for( size_t i = 0; i<Dimension; i++ ) {
		this->m_PointValues[i].set_size( this->m_NumberOfPoints );
		this->m_Phi.mult( coeff[i], this->m_PointValues[i] );
	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::InterpolateField() {
	const DimensionParameters coeff = this->VectorizeCoefficients();
	// Check m_Phi and initializations
	if( this->m_FieldPhi.rows() == 0 || this->m_FieldPhi.cols() == 0 ) {
		this->ComputeMatrix( Self::PHI_FIELD );
	}

	size_t npix = this->m_DisplacementField->GetLargestPossibleRegion().GetNumberOfPixels();

	DimensionVector interpField[Dimension];
	for( size_t i = 0; i<Dimension; i++ ) {
		interpField[i] = DimensionVector();
		interpField[i].set_size(npix);
		this->m_FieldPhi.mult( coeff[i], interpField[i] );
	}

	bool setVector;
	ScalarType val;
	VectorType v; v.Fill(0.0);

	FieldPointer field = FieldType::New();
	field->SetRegions( this->m_DisplacementField->GetLargestPossibleRegion().GetSize() );
	field->SetOrigin( this->m_DisplacementField->GetOrigin() );
	field->SetSpacing( this->m_DisplacementField->GetSpacing() );
	field->SetDirection( this->m_DisplacementField->GetDirection() );
	field->Allocate();
	field->FillBuffer( v );

	VectorType* obuf = field->GetBufferPointer();

	for( size_t row = 0; row<npix; row++ ) {
		v.Fill( 0.0 );
		setVector = false;
		for( size_t i = 0; i<Dimension; i++) {
			val = interpField[i][row];

			if( fabs(val) > 1.0e-5 ){
				v[i] = val;
				setVector = true;
			}
		}
		if (setVector) {
			*( obuf + row ) = v;
		}
	}
	this->SetDisplacementField( field );
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::ComputeInverse() {
	InvertFieldPointer invfield = InvertFieldFilter::New();
	invfield->SetDisplacementField(this->GetDisplacementField());
	invfield->Update();

	this->SetInverseDisplacementField(invfield->GetOutput());

	//// Check m_Phi_inverse and initializations
	//if( this->m_Phi_inverse.rows() == 0 || this->m_Phi_inverse.cols() == 0 ) {
	//	this->InvertPhi();
	//}
    //
	//DimensionParameters coeff;
	//// Compute coefficients
	//for( size_t i = 0; i<Dimension; i++ ) {
	//	coeff[i] = DimensionVector( this->m_NumberOfPoints );
	//	this->m_Phi_inverse.mult( this->m_PointValues[i], coeff[i] );
	//}
    //
	//// TODO set coeff here or inside Interpolate( coeff )
    //
	//// Interpolation with new coefficients
	//this->UpdateField( coeff );
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::UpdateField( const DimensionParameters& coeff ) {

	// TODO: this needs a re-implementation
	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeMatrix( Self::S );
	}

	// Interpolate
	DimensionParameters fieldValues;
	for( size_t i = 0; i<Dimension; i++)
		this->m_S.mult( coeff[i], fieldValues[i] );
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::ComputeCoefficients() {
	if ( this->m_NumberOfDimParameters == 0 ) {
		this->InitializeCoefficientsImages();
	}

	if( this->m_S.rows() == 0 || this->m_S.cols() == 0 ) {
		this->ComputeMatrix( Self::S );
	}

	SolverVector X[Dimension], Y[Dimension];

	DimensionParameters fieldValues = this->VectorizeField( this->m_DisplacementField );
	DimensionParameters coeffs;

	for ( size_t i = 0; i<Dimension; i++) {
		// TODO: check m_NumberOfPoints here
		// Y[i] = SolverVector( this->m_NumberOfPoints );
		Y[i] = SolverVector( fieldValues[i].size() );
		vnl_copy< DimensionVector, SolverVector >( fieldValues[i], Y[i] );
		X[i] = SolverVector( this->m_NumberOfDimParameters );
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
		size_t offset = col * this->m_NumberOfDimParameters;
		for( size_t k = 0; k<this->m_NumberOfDimParameters; k++) {
			this->m_Parameters[k + offset] = X[col][k];
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

	DimensionParameters coeff = this->VectorizeCoefficients();
	DimensionParameters result[Dimension];
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

	for ( size_t row = 0; row<this->m_NumberOfDimParameters; row++ ) {
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
const typename SparseMatrixTransform<TScalar,NDimensions>::WeightsMatrix*
SparseMatrixTransform<TScalar,NDimensions>
::GetPhi(const bool onlyvalid) {
	// Check m_Phi and initializations
	if( this->m_Phi.rows()==0 || this->m_Phi.cols()==0 ) {
		this->ComputeMatrix( Self::PHI );
	}

	if (onlyvalid) {
		return &this->m_Phi_valid;
	} else {
		return &this->m_Phi;
	}
}

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

	if (this->m_NumberOfDimParameters == 0) {
		this->m_ControlGridDirection = images[0]->GetDirection();
		this->m_ControlGridOrigin = images[0]->GetOrigin();
		this->m_ControlGridSize = images[0]->GetLargestPossibleRegion().GetSize();
		this->m_ControlGridSpacing = images[0]->GetSpacing();
		this->InitializeCoefficientsImages();
	}


	for( size_t dim = 0; dim < Dimension; dim++ ) {
		CoeffImageConstPointer c = itkDynamicCastInDebugMode< const CoefficientsImageType* >( images[dim].GetPointer() );
		this->SetCoefficientsImage( dim, c );
	}
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::SetCoefficientsVectorImage( const FieldType* images ) {
	this->SetCoefficientsField(images);
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::AddCoefficientsVectorImage( const FieldType* images ) {
	ScalarType* buff[Dimension];

	if (this->m_NumberOfDimParameters == 0) {
		this->SetCoefficientsField(images);
	} else {
		this->ComputeCoefficients();

		buff = this->m_CoefficientsField->GetBufferPointer();

		const VectorType* fbuf = images->GetBufferPointer();
		VectorType v;
		for( size_t row = 0; row < this->m_NumberOfDimParameters; row++ ) {
			*( buff + row )+= *( fbuf + row );
		}
	}
}


template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::SetCoefficientsImage( size_t dim, const CoefficientsImageType* c ) {
	if ( dim >= Dimension || dim < 0 ) {
		itkExceptionMacro(<< "trying to set coefficients for undefined dimension");
	}
	size_t offset = dim * this->m_NumberOfDimParameters;
	const ScalarType* fbuf = c->GetBufferPointer();
	for( size_t row = 0; row < this->m_NumberOfDimParameters; row++ ) {
		this->m_Parameters[row + offset] = *(fbuf + row);
	}
}


template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionParameters
SparseMatrixTransform<TScalar,NDimensions>
::VectorizeCoefficients() const {
	DimensionParameters m;
	const ScalarType* cbuffer;
	for( size_t i = 0; i<Dimension; i++) {
		m[i] = DimensionVector( this->m_NumberOfDimParameters );

		for( size_t j = 0; j < this->m_NumberOfDimParameters; j++ ){
			m[i][j] = this->m_Parameters[i * this->m_NumberOfDimParameters + j];
		}
	}
	return m;
}

template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionParameters
SparseMatrixTransform<TScalar,NDimensions>
::VectorizeDerivatives() const {
	DimensionParameters m;
	const ScalarType* cbuffer;
	for( size_t i = 0; i<Dimension; i++) {
		m[i] = DimensionVector( this->m_NumberOfDimParameters );
		cbuffer = this->m_Derivatives[i]->GetBufferPointer();

		for( size_t j = 0; j < this->m_NumberOfDimParameters; j++ ){
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

	for( size_t row = 0; row<this->m_NumberOfDimParameters; row++ ) {
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
	for( size_t row = 0; row<this->m_NumberOfDimParameters; row++ ) {
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
typename SparseMatrixTransform<TScalar,NDimensions>::DimensionParameters
SparseMatrixTransform<TScalar,NDimensions>
::VectorizeField( const FieldType* image ) {
	DimensionParameters vectorized;

	for( size_t col = 0; col<Dimension; col++) {
		vectorized[col] = DimensionVector( image->GetLargestPossibleRegion().GetNumberOfPixels() );
		vectorized[col].fill(0.0);
	}

	VectorType v;
	const VectorType *cbuffer = image->GetBufferPointer();

	for( size_t row = 0; row<this->m_NumberOfDimParameters; row++ ) {
		v = *( cbuffer + row );
		for( size_t col = 0; col<Dimension; col++) {
			vectorized[col][row] = v[col];
		}
	}
	return vectorized;
}

template< class TScalar, unsigned int NDimensions >
typename SparseMatrixTransform<TScalar,NDimensions>::AltCoeffPointer
SparseMatrixTransform<TScalar,NDimensions>
::GetFlatParameters() {
	this->m_FlatCoeffs = AltCoeffType::New();
	AltCoeffContainerPointer points = this->m_FlatCoeffs->GetPoints();
	AltCoeffDataPointer data = this->m_FlatCoeffs->GetPointData();

	for(size_t i = 0; i < this->m_NumberOfDimParameters; i++) {
		points->InsertElement(i, this->m_ParamLocations[i]);
		data->InsertElement(i, this->GetCoefficient(i));
	}
	return this->m_FlatCoeffs;
}

template< class TScalar, unsigned int NDimensions >
void
SparseMatrixTransform<TScalar,NDimensions>
::Initialize() {
	// Hypothesis 1: coefficients image reference has been set (origin, size, spacing)
	// Nothing to do

	// Hypothesis 2: m_DomainExtent and m_ControlGridSpacing are set
	double extent[Dimension];
	double center[Dimension];
	bool extentDefined = true;
	bool spacingDefined = true;
	for( size_t i = 0; i < Dimension; i++) {
		center[i] = 0.5 * (this->m_DomainExtent[1][i] + this->m_DomainExtent[0][i]);
		extent[i] = fabs(this->m_DomainExtent[1][i] - this->m_DomainExtent[0][i]);

		if (extent[i] < 1.0) {
			extentDefined = false;
		}

		if( fabs(this->m_ControlGridSpacing[i]) < 1.0) {
			spacingDefined = false;
		}
	}

	bool outputDefined = this->m_PointLocations.size() > 0 || this->m_ParamLocations.size() > 0;

	if(extentDefined && spacingDefined) {
		for(size_t i=0; i<Dimension; i++) {
			size_t half = floor(0.5*extent[i]/this->m_ControlGridSpacing[i]) + 1;
			this->m_ControlGridSize[i] = half * 2 + 1;
			this->m_ControlGridOrigin[i] = center[i] - (half * this->m_ControlGridSpacing[i]);
		}
	} else if (spacingDefined && outputDefined) {
		// Hypothesis 3: output (points or grid) and spacing have been set
		itkExceptionMacro(<< "not implemented");
	}

	this->InitializeCoefficientsImages();
};

}


#endif /* SPARSEMATRIXTRANSFORM_HXX_ */
