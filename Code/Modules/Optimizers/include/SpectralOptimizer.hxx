// --------------------------------------------------------------------------------------
// File:          SpectralOptimizer.hxx
// Date:          Jul 31, 2013
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

#ifndef SPECTRALOPTIMIZER_HXX_
#define SPECTRALOPTIMIZER_HXX_

#include "SpectralOptimizer.h"

#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <cmath>
#include <itkComplexToRealImageFilter.h>
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;

#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

namespace rstk {

/**
 * Default constructor
 */
template< typename TFunctional >
SpectralOptimizer<TFunctional>::SpectralOptimizer() {
	this->m_LearningRate = itk::NumericTraits<InternalComputationValueType>::One;
	this->m_CurrentValue = itk::NumericTraits<MeasureType>::infinity();
	this->m_MinimumConvergenceValue = 1e-8;
	this->m_ConvergenceWindowSize = 30;
	this->m_StepSize =  1.0;
	this->SetAlpha( 2.0 );
	this->SetBeta( 0.0 );

	this->m_NumberOfIterations = 100;
	this->m_CurrentIteration   = 0;
	this->m_StopCondition      = MAXIMUM_NUMBER_OF_ITERATIONS;
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_GridSize.Fill( 15 );
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
	os << indent << "Learning rate:" << this->m_LearningRate << std::endl;
	os << indent << "Number of iterations: " << this->m_NumberOfIterations << std::endl;
	os << indent << "Current iteration: " << this->m_CurrentIteration << std::endl;
	os << indent << "Stop condition:" << this->m_StopCondition << std::endl;
	os << indent << "Stop condition description: " << this->m_StopConditionDescription.str() << std::endl;
}



template< typename TFunctional >
void SpectralOptimizer<TFunctional>::Start() {
	itkDebugMacro("SpectralOptimizer::Start()");


	/* Settings validation */
	if ( this->m_Functional.IsNull() ) {
		itkExceptionMacro("Energy functional must be set");
	}

	/* Check & initialize parameter fields */
	this->InitializeParameters();

	this->InitializeAuxiliarParameters();

	/* Initialize convergence checker */
	this->m_ConvergenceMonitoring = ConvergenceMonitoringType::New();
	this->m_ConvergenceMonitoring->SetWindowSize( this->m_ConvergenceWindowSize );

//	if( this->m_ReturnBestParametersAndValue )	{
//		this->m_BestParameters = this->GetCurrentPosition( );
//		this->m_CurrentBestValue = NumericTraits< MeasureType >::max();
//	}

	this->m_Functional->CopyInformation( this->m_Parameters );
	this->m_Functional->Initialize();

	this->m_CurrentIteration = 0;
	this->Resume();
}

template< typename TFunctional >
const typename SpectralOptimizer<TFunctional>::StopConditionReturnStringType
SpectralOptimizer<TFunctional>::
GetStopConditionDescription() const {
  return this->m_StopConditionDescription.str();
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::Stop() {
	itkDebugMacro( "Stop called with a description - "
	  << this->GetStopConditionDescription() );
	this->m_Stop = true;
	this->InvokeEvent( itk::EndEvent() );

//	if( this->m_ReturnBestParametersAndValue )	{
//		this->GetMetric()->SetParameters( this->m_BestParameters );
//		this->m_CurrentValue = this->m_CurrentBestValue;
//	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::Resume() {
	this->m_StopConditionDescription.str("");
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->InvokeEvent( itk::StartEvent() );
	this->m_Stop = false;


	while( ! this->m_Stop )	{
		/* Compute functional value/derivative. */
		try	{
			this->m_Functional->ComputeDerivative();
		}
		catch ( itk::ExceptionObject & err ) {
			this->m_StopCondition = COSTFUNCTION_ERROR;
			this->m_StopConditionDescription << "Functional error during optimization";
			this->Stop();
			throw err;  // Pass exception to caller
		}

		/* Check if optimization has been stopped externally.
		 * (Presumably this could happen from a multi-threaded client app?) */
		if ( this->m_Stop ) {
			this->m_StopConditionDescription << "Stop() called externally";
			break;
		}


		/* Advance one step along the gradient.
		 * This will modify the gradient and update the transform. */
		this->Iterate();

		this->InvokeEvent( itk::IterationEvent() );


		/* Update the level sets contour and deformation field */
		this->m_Functional->UpdateContour( this->m_NextParameters );

		this->ComputeIterationChange();

		this->SetUpdate();


		/* TODO Store best value and position */
		//if ( this->m_ReturnBestParametersAndValue && this->m_CurrentValue < this->m_CurrentBestValue )
		//{
		//	this->m_CurrentBestValue = this->m_CurrentValue;
		//	this->m_BestParameters = this->GetCurrentPosition( );
		//}

#ifndef NDEBUG
		size_t nContours =this->m_Functional->GetCurrentContours().size();

		typedef rstk::DisplacementFieldComponentsFileWriter<ParametersType> Writer;
		typename Writer::Pointer p = Writer::New();
		std::stringstream ss;
		ss.str("");
		ss << "parameters_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
		p->SetFileName( ss.str().c_str() );
		p->SetInput( this->m_Parameters );
		p->Update();

		typedef itk::ImageFileWriter< typename FunctionalType::ProbabilityMapType > MapWriter;
		for( size_t r = 0; r <= nContours; r++){
			ss.str("");
			ss << "region_" << r << "_it" << std::setfill('0')<<std::setw(3) << this->m_CurrentIteration << ".nii.gz";
			typename MapWriter::Pointer wr = MapWriter::New();
			wr->SetInput( this->m_Functional->GetCurrentMap(r));
			wr->SetFileName(ss.str().c_str() );
			wr->Update();
		}
#endif


		std::cout << "[" << this->m_CurrentIteration << "] " << this->m_CurrentValue << " | " << this->m_Functional->GetValue()
				<< " | " << this->GetCurrentRegularizationEnergy() << std::endl;

		/*
		 * Check the convergence by WindowConvergenceMonitoringFunction.
		 */

		/*
		this->m_ConvergenceMonitoring->AddEnergyValue( this->m_CurrentValue );
		try {
			this->m_ConvergenceValue = this->m_ConvergenceMonitoring->GetConvergenceValue();
			if (this->m_ConvergenceValue <= this->m_MinimumConvergenceValue) {
				this->m_StopConditionDescription << "Convergence checker passed at iteration " << this->m_CurrentIteration << ".";
				this->m_StopCondition = Superclass::CONVERGENCE_CHECKER_PASSED;
				this->Stop();
				break;
			}
		}
		catch(std::exception & e) {
			std::cerr << "GetConvergenceValue() failed with exception: " << e.what() << std::endl;
		}


		if( this->m_CurrentValue < this->m_ConvergenceValue ) {
			this->m_StopConditionDescription << "Parameters field changed below the minimum threshold.";
			this->m_StopCondition = Superclass::STEP_TOO_SMALL;
			this->Stop();
			break;
		}
		*/

		/* Update and check iteration count */
		this->m_CurrentIteration++;

		if ( this->m_CurrentIteration >= this->m_NumberOfIterations ) {
			this->m_StopConditionDescription << "Maximum number of iterations (" << this->m_NumberOfIterations << ") exceeded.";
			this->m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
			this->Stop();
			break;
		}
	} //while (!m_Stop)
}


template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::GetCurrentRegularizationEnergy() {
	this->m_RegularizationEnergy=0;
	const VectorType* fBuffer = this->m_LastField->GetBufferPointer();
	size_t nPix = this->m_LastField->GetLargestPossibleRegion().GetNumberOfPixels();

	VectorType u;

	for ( size_t pix = 0; pix<nPix; pix++) {
		u = *(fBuffer+pix);
		this->m_RegularizationEnergy+= dot_product(u.GetVnlVector(), this->m_A.GetVnlMatrix() * u.GetVnlVector() );
	}

	double normalizer = 1.0;

	for( size_t i = 0; i<Dimension; i++ ) {
		normalizer*= this->m_LastField->GetSpacing()[i];
	}

	this->m_Functional->GetFieldInterpolator()->ComputeJacobian();

	typedef typename FunctionalType::FieldInterpolatorType::JacobianType JacobianType;
	JacobianType j;
	JacobianType M;

	for ( size_t pix = 0; pix<nPix; pix++) {
		j = this->m_Functional->GetFieldInterpolator()->GetJacobian(pix);
		M = (j * m_B) * j.GetTranspose();
		for ( size_t i = 0; i<Dimension; i++ )
			this->m_RegularizationEnergy+= M(i,i);
	}


	return normalizer * this->m_RegularizationEnergy;
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::SpectralUpdate(
		ParametersPointer parameters,
		const ParametersType* lambda,
		ParametersPointer nextParameters,
		bool changeDirection )
{
	itkDebugMacro("Optimizer Spectral Update");
	InternalComputationValueType min_val = 1.0e-5;
	InternalComputationValueType dir = (changeDirection)?(-1.0):1.0;

	// Get deformation fields pointers
	const VectorType* uFieldBuffer = parameters->GetBufferPointer();
	const VectorType* lFieldBuffer = lambda->GetBufferPointer();
	VectorType* nextBuffer = nextParameters->GetBufferPointer();

	// Create container for deformation field components
	ParametersComponentPointer fieldComponent = ParametersComponentType::New();
	fieldComponent->SetDirection( this->m_Parameters->GetDirection() );
	fieldComponent->SetRegions( this->m_Parameters->GetLargestPossibleRegion() );
	fieldComponent->SetSpacing( this->m_Parameters->GetSpacing() );
	fieldComponent->SetOrigin( this->m_Parameters->GetOrigin());
	fieldComponent->Allocate();
	PointValueType* dCompBuffer = fieldComponent->GetBufferPointer();

	ComplexFieldPointer fftNum;

	size_t nPix = fieldComponent->GetLargestPossibleRegion().GetNumberOfPixels();
	InternalComputationValueType r = 1.0/this->m_StepSize;
	VectorType u,l;

	for( size_t d = 0; d < Dimension; d++ ) {

		for(size_t pix = 0; pix< nPix; pix++) {
			u = *(uFieldBuffer+pix);
			l = *(lFieldBuffer+pix);
			*(dCompBuffer+pix) =  u[d] * r + dir * l[d];
		}

		FFTPointer fftFilter = FFTType::New();
		fftFilter->SetInput( fieldComponent );
		fftFilter->Update();  // This is required for computing the denominator for first time
		FTDomainPointer fftComponent =  fftFilter->GetOutput();

		if ( fftNum.IsNull() ) {
			fftNum = ComplexFieldType::New();
			fftNum->SetRegions( fftComponent->GetLargestPossibleRegion() );
			fftNum->Allocate();
			fftNum->FillBuffer( ComplexType( itk::NumericTraits<ComplexValueType>::Zero, itk::NumericTraits<ComplexValueType>::Zero) );
		}

		// Fill in ND fourier transform of numerator
		ComplexFieldValue* fftBuffer = fftNum->GetBufferPointer();
		ComplexType* fftCompBuffer = fftComponent->GetBufferPointer();
		size_t nFrequel = fftComponent->GetLargestPossibleRegion().GetNumberOfPixels();
		ComplexFieldValue cur;
		for(size_t fr = 0; fr< nFrequel; fr++) {
			cur = *(fftBuffer+fr);
			cur[d] = *(fftCompBuffer+fr);
			*(fftBuffer+fr) = cur;
		}
	}

	this->ApplyRegularizationTerm( fftNum );

	for( size_t d = 0; d < Dimension; d++ ) {
		// Take out fourier component
		FTDomainPointer fftComponent = FTDomainType::New();
		fftComponent->SetRegions( fftNum->GetLargestPossibleRegion() );
		fftComponent->SetOrigin( fftNum->GetOrigin() );
		fftComponent->SetSpacing( fftNum->GetSpacing() );
		fftComponent->Allocate();
		fftComponent->FillBuffer( ComplexType(0.0,0.0) );
		ComplexType* fftCompBuffer = fftComponent->GetBufferPointer();
		ComplexFieldValue* fftBuffer = fftNum->GetBufferPointer();
		size_t nFrequel = fftComponent->GetLargestPossibleRegion().GetNumberOfPixels();
		for(size_t fr = 0; fr< nFrequel; fr++) {
			*(fftCompBuffer+fr) = (*(fftBuffer+fr))[d];
		}

		IFFTPointer ifft = IFFTType::New();
		ifft->SetInput( fftComponent );

		try {
			ifft->Update();
		}
		catch ( itk::ExceptionObject & err ) {
			this->m_StopCondition = UPDATE_PARAMETERS_ERROR;
			this->m_StopConditionDescription << "Optimizer error: spectral update failed";
			this->Stop();
			// Pass exception to caller
			throw err;
		}

		// Set component on this->m_ParametersNext
		const PointValueType* resultBuffer = ifft->GetOutput()->GetBufferPointer();
		PointValueType val;
		for(size_t pix = 0; pix< nPix; pix++) {
			val = *(resultBuffer+pix);
			(*(nextBuffer+pix))[d] = ( fabs(val) > min_val )?val:0.0;
		}
	}

}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::ApplyRegularizationTerm( typename SpectralOptimizer<TFunctional>::ComplexFieldType* reference ){
	if( this->m_Denominator.IsNull() ) { // If executed for first time, cache the map
		this->InitializeDenominator( reference );
	}

	// Fill reference buffer with elements updated with the filter
	ComplexFieldValue* nBuffer = reference->GetBufferPointer();
	MatrixType* dBuffer = this->m_Denominator->GetBufferPointer();
	size_t nPix = this->m_Denominator->GetLargestPossibleRegion().GetNumberOfPixels();
	ComplexFieldValue curval;
	VectorType realval;
	VectorType imagval;
	MatrixType regTensor;

	for (size_t pix = 0; pix < nPix; pix++ ) {
		curval = *(nBuffer+pix);
		regTensor = *(dBuffer+pix);

		for( size_t i = 0; i<Dimension; i++ ) {
			realval[i] = curval[i].real();
			imagval[i] = curval[i].imag();
		}

		realval = regTensor * realval;
		imagval = regTensor * imagval;

		for( size_t i = 0; i<Dimension; i++ ) {
			curval[i] = ComplexType(realval[i], imagval[i]);
		}


		*(nBuffer+pix) = curval;
	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::InitializeDenominator( typename SpectralOptimizer<TFunctional>::ComplexFieldType* reference ){
	InternalComputationValueType pi2 = 2* vnl_math::pi;
	InternalComputationValueType tfFactor = this->m_Parameters->GetLargestPossibleRegion().GetNumberOfPixels();
	InternalComputationValueType singular_threshold = 1.0e-9;
	MatrixType I, t;
	I.SetIdentity();

	this->m_Denominator = TensorFieldType::New();
	this->m_Denominator->SetSpacing(   reference->GetSpacing() );
	this->m_Denominator->SetDirection( reference->GetDirection() );
	this->m_Denominator->SetRegions( reference->GetLargestPossibleRegion().GetSize() );
	this->m_Denominator->Allocate();

	MatrixType* buffer = this->m_Denominator->GetBufferPointer();

	MatrixType C = ( I * (1.0/this->m_StepSize) + m_A ) * tfFactor;

	// Check that C is not singular
	vnl_matrix_inverse< ComplexValueType > inv_C( C.GetVnlMatrix() );
	InternalComputationValueType detC = inv_C.determinant_magnitude();

	if( detC <= singular_threshold ) {
		itkExceptionMacro( << "norm penalizer is singular.")
	}

	t = MatrixType(inv_C);
	this->m_Denominator->FillBuffer( t );

	size_t nPix = this->m_Denominator->GetLargestPossibleRegion().GetNumberOfPixels();

	// Check that B is not singular
	vnl_matrix_inverse< ComplexValueType > inv_B( m_B.GetVnlMatrix() );
	InternalComputationValueType detB = inv_B.determinant_magnitude();

	if( detB > singular_threshold ) {
		typename TensorFieldType::SizeType size = this->m_Denominator->GetLargestPossibleRegion().GetSize();

		// Fill Buffer with denominator data
		InternalComputationValueType lag_el; // accumulates the FT{lagrange operator}.
		typename TensorFieldType::IndexType idx;
		for (size_t pix = 0; pix < nPix; pix++ ) {
			lag_el = 0.0;
			idx = this->m_Denominator->ComputeIndex( pix );
			for(size_t d = 0; d < Dimension; d++ ) {
				lag_el+= 2.0 * cos( (pi2*idx[d])/size[d]) - 2.0;
			}
			t = C - m_B * lag_el * (-1.0) * tfFactor;

			*(buffer+pix) = t.GetInverse();
		}
	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::ComputeIterationChange() {
	const VectorType* fnextBuffer = this->m_Functional->GetCurrentDisplacementField()->GetBufferPointer();

#ifndef NDEBUG
		typedef rstk::DisplacementFieldFileWriter<typename FunctionalType::FieldType> Writer;
		typename Writer::Pointer p = Writer::New();
		std::stringstream ss;
		ss.str("");
		ss << "field_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
		p->SetFileName( ss.str().c_str() );
		p->SetInput( this->m_Functional->GetCurrentDisplacementField() );
		p->Update();
#endif

	VectorType* fBuffer = this->m_LastField->GetBufferPointer();
	size_t nPix = this->m_LastField->GetLargestPossibleRegion().GetNumberOfPixels();

	InternalComputationValueType totalNorm = 0;
	VectorType t0,t1;
	InternalComputationValueType diff = 0.0;
	for (size_t pix = 0; pix < nPix; pix++ ) {
		t0 = *(fBuffer+pix);
		t1 = *(fnextBuffer+pix);
		diff = ( t1 - t0 ).GetNorm();
		totalNorm += diff;

		*(fBuffer+pix) = t1; // Copy current to last, once evaluated
	}

	this->m_CurrentValue = totalNorm/nPix;
}


template< typename TFunctional >
void SpectralOptimizer<TFunctional>::InitializeParameters() {
	VectorType zerov = itk::NumericTraits<VectorType>::Zero;

	if (this->m_Parameters.IsNull()) {
		PointType orig = this->m_Functional->GetReferenceImage()->GetOrigin();
		PointType end;
		typename ParametersType::IndexType end_idx;
		typename ParametersType::SizeType refSize = this->m_Functional->GetReferenceImage()->GetLargestPossibleRegion().GetSize();

		for( size_t dim = 0; dim < Dimension; dim++ ) end_idx[dim] = refSize[dim];
		this->m_Functional->GetReferenceImage()->TransformIndexToPhysicalPoint( end_idx, end );

		this->m_Parameters = ParametersType::New();

		typename ParametersType::SpacingType spacing;
		for( size_t dim = 0; dim < Dimension; dim++ ) {
			spacing[dim] = fabs( (end[dim]-orig[dim])/(1.0*this->m_GridSize[dim]));
		}

		VectorType zerov;
		zerov.Fill(0.0);

		this->m_Parameters->SetRegions( this->m_GridSize );
		this->m_Parameters->SetSpacing( spacing );
		this->m_Parameters->SetDirection( this->m_Functional->GetReferenceImage()->GetDirection() );
		this->m_Parameters->SetOrigin( orig );
		this->m_Parameters->Allocate();
		this->m_Parameters->FillBuffer( zerov );

	}

	/* Initialize next parameters */
	this->m_NextParameters = ParametersType::New();
	this->m_NextParameters->SetRegions( this->m_Parameters->GetLargestPossibleRegion() );
	this->m_NextParameters->SetSpacing( this->m_Parameters->GetSpacing() );
	this->m_NextParameters->SetDirection( this->m_Parameters->GetDirection() );
	this->m_NextParameters->SetOrigin( this->m_Parameters->GetOrigin() );
	this->m_NextParameters->Allocate();
	this->m_NextParameters->FillBuffer( zerov );


	this->m_LastField = ParametersType::New();
	this->m_LastField->SetRegions( this->m_Parameters->GetLargestPossibleRegion() );
	this->m_LastField->SetSpacing( this->m_Parameters->GetSpacing() );
	this->m_LastField->SetDirection( this->m_Parameters->GetDirection() );
	this->m_LastField->SetOrigin( this->m_Parameters->GetOrigin() );
	this->m_LastField->Allocate();
	this->m_LastField->FillBuffer( zerov );
	//this->m_CurrentField = ParametersType::New();
	//this->m_CurrentField->SetRegions( this->m_Parameters->GetLargestPossibleRegion() );
	//this->m_CurrentField->SetSpacing( this->m_Parameters->GetSpacing() );
	//this->m_CurrentField->SetDirection( this->m_Parameters->GetDirection() );
	//this->m_CurrentField->SetOrigin( this->m_Parameters->GetOrigin() );
	//this->m_CurrentField->Allocate();
	//this->m_CurrentField->FillBuffer( zerov );
}

} // end namespace rstk


#endif /* SPECTRALOPTIMIZER_HXX_ */
