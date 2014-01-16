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
	this->m_ConvergenceWindowSize = 50;
	this->m_StepSize =  1.0;
	this->m_Alpha.Fill( 1.0 );
	this->m_Beta.Fill( 10.0 );

	this->m_NumberOfIterations = 250;
	this->m_CurrentIteration   = 0;
	this->m_StopCondition      = MAXIMUM_NUMBER_OF_ITERATIONS;
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_GridSize.Fill( 20 );

	m_CurrentTotalEnergy = itk::NumericTraits<MeasureType>::infinity();
	m_RegularizationEnergyUpdated = false;
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
	this->InvokeEvent( itk::StartEvent() );
	this->m_CurrentIteration = 1;
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
		this->m_CurrentSpeeds = this->m_Functional->GetDerivative();

#ifndef NDEBUG
		typedef rstk::DisplacementFieldComponentsFileWriter<ParametersType> ComponentsWriter;
		typename ComponentsWriter::Pointer p = ComponentsWriter::New();
		std::stringstream ss;
		ss.str("");
		ss << "speeds_" << std::setfill('0')  << std::setw(3) << this->GetCurrentIteration();
		p->SetFileName( ss.str().c_str() );
		p->SetInput( this->m_CurrentSpeeds );
		p->Update();
#endif


		this->Iterate();

		/* Update the deformation field */
		this->UpdateField();
		this->m_CurrentValue = this->ComputeIterationChange();
		this->SetUpdate();


		/* TODO Store best value and position */
		//if ( this->m_ReturnBestParametersAndValue && this->m_CurrentValue < this->m_CurrentBestValue )
		//{
		//	this->m_CurrentBestValue = this->m_CurrentValue;
		//	this->m_BestParameters = this->GetCurrentPosition( );
		//}


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
		if ( this->m_CurrentIteration >= this->m_NumberOfIterations ) {
			this->m_StopConditionDescription << "Maximum number of iterations (" << this->m_NumberOfIterations << ") exceeded.";
			this->m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
			this->Stop();
			break;
		}


		this->InvokeEvent( itk::IterationEvent() );
		this->m_CurrentIteration++;
	} //while (!m_Stop)
}


template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::GetCurrentRegularizationEnergy() {
	if (!this->m_RegularizationEnergyUpdated ){
		//this->m_RegularizationEnergy=0;
		//const VectorType* fBuffer = this->m_LastField->GetBufferPointer();
		//size_t nPix = this->m_LastField->GetLargestPossibleRegion().GetNumberOfPixels();
        //
		//VectorType u;
        //
		//for ( size_t pix = 0; pix<nPix; pix++) {
		//	u = *(fBuffer+pix);
		//	for ( size_t i = 0; i<Dimension; i++) {
		//		u[i] = u[i]*u[i];
		//	}
		//	this->m_RegularizationEnergy+= this->m_Alpha * u;
		//}
        //
		//double normalizer = 1.0;
        //
		//for( size_t i = 0; i<Dimension; i++ ) {
		//	normalizer*= this->m_LastField->GetSpacing()[i];
		//}
        //
		//this->m_Functional->GetFieldInterpolator()->ComputeJacobian();
        //
		//typedef typename FunctionalType::FieldInterpolatorType::JacobianType JacobianType;
		//JacobianType j;
        //
		//for ( size_t pix = 0; pix<nPix; pix++) {
		//	j = this->m_Functional->GetFieldInterpolator()->GetJacobian(pix);
		//	for ( size_t i = 0; i<Dimension; i++) {
		//		u[i] = j[i][i]*j[i][i];
		//	}
		//	this->m_RegularizationEnergy+= this->m_Beta * u;
		//}
        //
		//this->m_RegularizationEnergy = normalizer * this->m_RegularizationEnergy;
		//this->m_RegularizationEnergyUpdated = true;
	}
	return this->m_RegularizationEnergy;
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::SpectralUpdate(
		ParametersType* parameters,
		const ParametersType* lambda,
		ParametersType* nextParameters,
		bool changeDirection )
{
	itkDebugMacro("Optimizer Spectral Update");

	double min_val = 1.0e-9;
	double dir = (changeDirection)?-1.0:1.0;
	double r = 1.0/this->m_StepSize;

	size_t nPix = parameters->GetLargestPossibleRegion().GetNumberOfPixels();


	std::vector< ParametersComponentPointer > numComps;

	for( size_t d = 0; d<Dimension; d++ ) {
		// Create container for deformation field components
		ParametersComponentPointer num = ParametersComponentType::New();
		num->SetDirection( this->m_Parameters->GetDirection() );
		num->SetRegions( this->m_Parameters->GetLargestPossibleRegion() );
		num->SetSpacing( this->m_Parameters->GetSpacing() );
		num->SetOrigin( this->m_Parameters->GetOrigin());
		num->Allocate();
		num->FillBuffer(0.0);
		numComps.push_back( num );
	}

	// Prepare numerator
	// Get deformation fields pointers
	const VectorType* uFieldBuffer = parameters->GetBufferPointer();
	const VectorType* lFieldBuffer = lambda->GetBufferPointer();
	VectorType u, l, n;

	for(size_t pix = 0; pix< nPix; pix++) {
		u = *(uFieldBuffer+pix);
		l = *(lFieldBuffer+pix);
		n = u * r + dir * l;

		if ( n.GetNorm()> 0.0 ) {

			for( size_t d = 0; d<Dimension; d++ ) {
				PointValueType* compBuffer = numComps[d]->GetBufferPointer();
				*(compBuffer+pix) = n[d];
			}
		}
	}

	VectorType* nextBuffer = nextParameters->GetBufferPointer();

	for( size_t d = 0; d < Dimension; d++ ) {
		FFTPointer fftFilter = FFTType::New();
		fftFilter->SetInput( numComps[d] );
		fftFilter->Update();  // This is required for computing the denominator for first time
		FTDomainPointer fftComponent =  fftFilter->GetOutput();

		this->ApplyRegularizationComponent(d, fftComponent );

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
::ApplyRegularizationComponent( size_t d, typename SpectralOptimizer<TFunctional>::FTDomainType* reference ){

	if( this->m_Denominator.size() != Dimension ) { // If executed for first time, cache the map
		this->InitializeDenominator( reference );
	}

	size_t nPix = this->m_Denominator[d]->GetLargestPossibleRegion().GetNumberOfPixels();

	// Fill reference buffer with elements updated with the filter
	ComplexType* nBuffer = reference->GetBufferPointer();
	PointValueType* dBuffer = this->m_Denominator[d]->GetBufferPointer();

	ComplexType curval;
	ComplexType ddor;
	ComplexType res;

	for (size_t pix = 0; pix < nPix; pix++ ) {
		curval = *(nBuffer+pix);
		ddor = ComplexType( *(dBuffer+pix), 0.0 );
		res = curval * ddor;
		*(nBuffer+pix) = curval * ddor;
	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::InitializeDenominator( itk::ImageBase<Dimension>* reference ){
	PointValueType pi2 = 2* vnl_math::pi;
	size_t nPix = reference->GetLargestPossibleRegion().GetNumberOfPixels();
	typename ParametersComponentType::SizeType size = reference->GetLargestPossibleRegion().GetSize();

	for (size_t i = 0; i<Dimension; i++) {
		PointValueType initVal = (1.0/this->m_StepSize) + this->m_Alpha[i];

		ParametersComponentPointer dimDdor = ParametersComponentType::New();
		dimDdor->SetSpacing(   reference->GetSpacing() );
		dimDdor->SetDirection( reference->GetDirection() );
		dimDdor->SetRegions( reference->GetLargestPossibleRegion().GetSize() );
		dimDdor->Allocate();

		PointValueType* buffer = dimDdor->GetBufferPointer();

		if( this->m_Beta.GetNorm() > 0.0 ) {
			// Fill Buffer with denominator data
			PointValueType lag_el; // accumulates the FT{lagrange operator}.
			typename ParametersComponentType::IndexType idx;
			for (size_t pix = 0; pix < nPix; pix++ ) {
				lag_el = 0.0;
				idx = dimDdor->ComputeIndex( pix );
				for(size_t d = 0; d < Dimension; d++ ) {
					lag_el+= 2.0 * cos( (pi2*idx[d])/size[d] ) - 2.0;
				}
				*(buffer+pix) = 1.0 / (initVal - m_Beta[i] * lag_el);
			}
		}
		else {
			dimDdor->FillBuffer( 1.0 / initVal );
		}


		this->m_Denominator.push_back( dimDdor );
	}
}


template< typename TFunctional >
void
SpectralOptimizer<TFunctional>::UpdateField() {
	bool update = false;

#ifndef NDEBUG
	typedef rstk::DisplacementFieldComponentsFileWriter<ParametersType> ComponentsWriter;
	typename ComponentsWriter::Pointer p = ComponentsWriter::New();
	std::stringstream ss;
	ss.str("");
	ss << "coefficients_" << std::setfill('0')  << std::setw(3) << this->GetCurrentIteration();
	p->SetFileName( ss.str().c_str() );
	p->SetInput( this->m_NextParameters );
	p->Update();
#endif

	FieldInterpolatorPointer fintp = this->m_Functional->GetFieldInterpolator();
	fintp->SetCoefficientsVectorImage( this->m_NextParameters );

	if (update) {
		fintp->UpdateField();
		typedef typename FieldInterpolatorType::FieldType FieldType;
		typename FieldType::ConstPointer f = fintp->GetField();

		itk::ImageAlgorithm::Copy<FieldType, FieldType>(f, this->m_CurrentDisplacementField,
			f->GetLargestPossibleRegion(),
			this->m_CurrentDisplacementField->GetLargestPossibleRegion()
	    );
	}
}

template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::ComputeIterationChange() {
	const VectorType* fnextBuffer = this->m_CurrentDisplacementField->GetBufferPointer();
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

	this->m_RegularizationEnergyUpdated = (totalNorm==0);

	return totalNorm/nPix;
}


template< typename TFunctional >
void SpectralOptimizer<TFunctional>::InitializeParameters() {
	VectorType zerov;
	zerov.Fill(0.0);

	if (this->m_Parameters.IsNull()) {
		ContinuousIndexType o_idx;
		o_idx.Fill( -0.5 );

		ContinuousIndexType e_idx;
		for ( size_t dim=0; dim< Dimension; dim++ ) {
			e_idx[dim] = this->m_Functional->GetReferenceImage()->GetLargestPossibleRegion().GetSize()[dim] - 0.5;
		}

		PointType first;
		PointType last;

		this->m_Functional->GetReferenceImage()->TransformContinuousIndexToPhysicalPoint( o_idx, first );
		this->m_Functional->GetReferenceImage()->TransformContinuousIndexToPhysicalPoint( e_idx, last );

		this->m_Parameters = ParametersType::New();


		PointType orig,step;
		typename ParametersType::SpacingType spacing;
		for( size_t dim = 0; dim < Dimension; dim++ ) {
			step[dim] = (last[dim]-first[dim])/(1.0*this->m_GridSize[dim]);
			spacing[dim] = fabs(step[dim]);
			orig[dim] = first[dim] + 0.5 * step[dim];
		}

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

	this->m_CurrentDisplacementField = ParametersType::New();
	this->m_CurrentDisplacementField->SetRegions( this->m_Parameters->GetLargestPossibleRegion() );
	this->m_CurrentDisplacementField->SetSpacing( this->m_Parameters->GetSpacing() );
	this->m_CurrentDisplacementField->SetDirection( this->m_Parameters->GetDirection() );
	this->m_CurrentDisplacementField->SetOrigin( this->m_Parameters->GetOrigin() );
	this->m_CurrentDisplacementField->Allocate();
	this->m_CurrentDisplacementField->FillBuffer( zerov );
}

template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::GetCurrentEnergy() {
	this->m_CurrentTotalEnergy = this->m_Functional->GetValue() + this->GetCurrentRegularizationEnergy();
	return this->m_CurrentTotalEnergy;
}

} // end namespace rstk


#endif /* SPECTRALOPTIMIZER_HXX_ */
