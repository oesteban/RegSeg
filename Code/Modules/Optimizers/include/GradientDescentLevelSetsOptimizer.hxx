// --------------------------------------------------------------------------
// File:             GradientDescentLevelSetsOptimizer.hxx
// Date:             01/11/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWEReg
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef GRADIENTDESCENTLEVELSETSOPTIMIZER_HXX_
#define GRADIENTDESCENTLEVELSETSOPTIMIZER_HXX_


#include "GradientDescentLevelSetsOptimizer.h"
#include <vector>
#include <vnl/vnl_math.h>
#include <cmath>
#include <itkComplexToRealImageFilter.h>
#include "DisplacementFieldFileWriter.h"


namespace rstk {

/**
 * Default constructor
 */
template< typename TLevelSetsFunction >
GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::GradientDescentLevelSetsOptimizer() {
	this->m_LearningRate = itk::NumericTraits<InternalComputationValueType>::One;
	this->m_MaximumStepSizeInPhysicalUnits = itk::NumericTraits<InternalComputationValueType>::Zero;
	this->m_MinimumConvergenceValue = 1e-3;
	this->m_ConvergenceWindowSize = 4;
	this->m_StepSize = 0.1;
	this->m_Alpha = 1.0;
	this->m_Beta = 1.0;
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
	os << indent << "Learning rate:" << this->m_LearningRate << std::endl;
	os << indent << "MaximumStepSizeInPhysicalUnits: "
             << this->m_MaximumStepSizeInPhysicalUnits << std::endl;
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::Start() {
	itkDebugMacro("GradientDescentLevelSetsOptimizer::Start()");

	VectorType zerov = itk::NumericTraits<VectorType>::Zero;
	/* Settings validation */
	if ( this->m_LevelSetsFunction.IsNull() ) {
		itkExceptionMacro("LevelSets Object must be set");
	}

	if (this->m_DeformationField.IsNull()) {
		/* TODO Initialize default deformation field */
		itkExceptionMacro("DeformationField Object must be set");
	}

	if ( this->m_SpeedsField.IsNull() ) {
		this->m_SpeedsField = DeformationFieldType::New();
		this->m_SpeedsField->SetRegions( this->m_DeformationField->GetLargestPossibleRegion() );
		this->m_SpeedsField->SetSpacing( this->m_DeformationField->GetSpacing() );
		this->m_SpeedsField->SetDirection( this->m_DeformationField->GetDirection() );
		this->m_SpeedsField->Allocate();
		this->m_SpeedsField->FillBuffer( zerov );
	}

	/* Initialize next deformationfield */
	this->m_NextDeformationField = DeformationFieldType::New();
	this->m_NextDeformationField->SetRegions( this->m_DeformationField->GetLargestPossibleRegion() );
	this->m_NextDeformationField->SetSpacing( this->m_DeformationField->GetSpacing() );
	this->m_NextDeformationField->SetDirection( this->m_DeformationField->GetDirection() );
	this->m_NextDeformationField->Allocate();
	this->m_NextDeformationField->FillBuffer( zerov );

	/* Initialize convergence checker */
	this->m_ConvergenceMonitoring = ConvergenceMonitoringType::New();
	this->m_ConvergenceMonitoring->SetWindowSize( this->m_ConvergenceWindowSize );

	/* Must call the superclass version for basic validation and setup */
	Superclass::Start();

//	if( this->m_ReturnBestParametersAndValue )	{
//		this->m_BestParameters = this->GetCurrentPosition( );
//		this->m_CurrentBestValue = NumericTraits< MeasureType >::max();
//	}

	this->m_CurrentIteration = 0;
	this->Resume();
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::Stop() {

//	if( this->m_ReturnBestParametersAndValue )	{
//		this->GetMetric()->SetParameters( this->m_BestParameters );
//		this->m_CurrentLevelSetsValue = this->m_CurrentBestValue;
//	}
	Superclass::Stop();
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::Resume() {
	this->m_StopConditionDescription.str("");
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->InvokeEvent( itk::StartEvent() );

	this->m_Stop = false;
	while( ! this->m_Stop )	{
		/* Compute metric value/derivative. */
		try	{
			this->m_SpeedsField = this->m_LevelSetsFunction->GetLevelSetsMap(this->m_DeformationField);
			/*
			typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;
			typename Writer::Pointer p = Writer::New();
			std::stringstream ss;
			ss << "speedsfield" << this->m_CurrentIteration << ".nii.gz";
			p->SetFileName( ss.str().c_str() );
			p->SetInput( this->m_SpeedsField );
			p->Update();
			*/
		}
		catch ( itk::ExceptionObject & err ) {
			this->m_StopCondition = Superclass::COSTFUNCTION_ERROR;
			this->m_StopConditionDescription << "LevelSets error during optimization";
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

		this->ComputeIterationEnergy();

		/*
		 * Check the convergence by WindowConvergenceMonitoringFunction.
		 */
		this->m_ConvergenceMonitoring->AddEnergyValue( this->m_CurrentLevelSetsValue );
		try
		{
			this->m_ConvergenceValue = this->m_ConvergenceMonitoring->GetConvergenceValue();
			if (this->m_ConvergenceValue <= this->m_MinimumConvergenceValue)
			{
				this->m_StopConditionDescription << "Convergence checker passed at iteration " << this->m_CurrentIteration << ".";
				this->m_StopCondition = Superclass::CONVERGENCE_CHECKER_PASSED;
				this->Stop();
				break;
			}
		}
		catch(std::exception & e) {
			std::cerr << "GetConvergenceValue() failed with exception: " << e.what() << std::endl;
		}

		/* Update the level sets contour and deformation field */
		this->m_LevelSetsFunction->UpdateDeformationField( this->m_NextDeformationField );
		itk::ImageAlgorithm::Copy<DeformationFieldType,DeformationFieldType>(
				this->m_NextDeformationField, this->m_DeformationField,
				this->m_NextDeformationField->GetLargestPossibleRegion(),
				this->m_DeformationField->GetLargestPossibleRegion());


		/* TODO Store best value and position */
		//if ( this->m_ReturnBestParametersAndValue && this->m_CurrentLevelSetsValue < this->m_CurrentBestValue )
		//{
		//	this->m_CurrentBestValue = this->m_CurrentLevelSetsValue;
		//	this->m_BestParameters = this->GetCurrentPosition( );
		//}

		std::cout << "Iteration " << this->m_CurrentIteration << std::endl;

		/* Update and check iteration count */
		this->m_CurrentIteration++;

		if ( this->m_CurrentIteration >= this->m_NumberOfIterations ) {
			this->m_StopConditionDescription << "Maximum number of iterations (" << this->m_NumberOfIterations << ") exceeded.";
			this->m_StopCondition = Superclass::MAXIMUM_NUMBER_OF_ITERATIONS;
			this->Stop();
			break;
		}
	} //while (!m_Stop)
}


template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::Iterate() {
	itkDebugMacro("Optimizer Iteration");

	// Get deformation fields pointers
	const VectorType* dFieldBuffer = this->m_DeformationField->GetBufferPointer();
	const VectorType* sFieldBuffer = this->m_SpeedsField->GetBufferPointer();
	VectorType* nextFieldBuffer = this->m_NextDeformationField->GetBufferPointer();

	// Create container for deformation field components
	DeformationComponentPointer fieldComponent = DeformationComponentType::New();
	fieldComponent->SetRegions( this->m_DeformationField->GetLargestPossibleRegion() );
	fieldComponent->SetSpacing( this->m_DeformationField->GetSpacing() );
	fieldComponent->Allocate();
	PointValueType* dCompBuffer = fieldComponent->GetBufferPointer();

	size_t nPix = fieldComponent->GetLargestPossibleRegion().GetNumberOfPixels();
	for( size_t d = 0; d < Dimension; d++ ) {
		// Get component from vector image
		size_t comp, idx2;
		for(size_t pix = 0; pix< nPix; pix++) {
			*(dCompBuffer+pix) = (PointValueType) (*(dFieldBuffer+pix))[d] / this->m_StepSize - (*(sFieldBuffer+pix))[d];
		}

		FFTPointer fft = FFTType::New();
		fft->SetInput( fieldComponent );
		fft->Update();  // This is required for computing the denominator for first time
		FTDomainPointer num = fft->GetOutput();
		this->ApplyRegularizationTerm( num );

		IFFTPointer ifft = IFFTType::New();
		ifft->SetInput( num );

		try {
			ifft->Update();
		}
		catch ( itk::ExceptionObject & err ) {
			this->m_StopCondition = Superclass::UPDATE_PARAMETERS_ERROR;
			this->m_StopConditionDescription << "UpdateDeformationField error";
			this->Stop();
			// Pass exception to caller
			throw err;
		}

		// Set component on this->m_NextDeformationField
		const PointValueType* resultBuffer = ifft->GetOutput()->GetBufferPointer();
		for(size_t pix = 0; pix< nPix; pix++) {
			(*(nextFieldBuffer+pix))[d] = *(resultBuffer+pix);
		}
	}

	this->InvokeEvent( itk::IterationEvent() );
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>
::ApplyRegularizationTerm( typename GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::FTDomainType* reference ){
	double pi2 = 2* vnl_math::pi;
	ComplexValueType constant = (1.0/this->m_StepSize) + m_Alpha;

	if( this->m_Denominator.IsNull() ) { // If executed for first time, cache the map
		this->m_Denominator = RealPartType::New();
		this->m_Denominator->SetRegions( reference->GetLargestPossibleRegion() );
		this->m_Denominator->SetSpacing( reference->GetSpacing() );
		this->m_Denominator->SetDirection( reference->GetDirection() );
		this->m_Denominator->Allocate();
		this->m_Denominator->FillBuffer( itk::NumericTraits<ComplexValueType>::Zero );

		// Fill Buffer with denominator data
		ComplexValueType lag_el; // accumulates the FT{lagrange operator}.
		typename RealPartType::IndexType idx;
		typename RealPartType::SizeType size = this->m_Denominator->GetLargestPossibleRegion().GetSize();
		ComplexValueType* buffer = this->m_Denominator->GetBufferPointer();
		size_t nPix = this->m_Denominator->GetLargestPossibleRegion().GetNumberOfPixels();
		for (size_t pix = 0; pix < nPix; pix++ ) {
			lag_el = 0.0;
			idx = this->m_Denominator->ComputeIndex( pix );
			for(size_t d = 0; d < Dimension; d++ ) lag_el+= 2.0*cos( pi2* idx[d]/(1.0*size[d]))-2;
			*(buffer+pix) = 1.0 / (constant - m_Beta* lag_el);
		}
	}

	// Fill reference buffer with elements updated with the filter
	ComplexType* nBuffer = reference->GetBufferPointer();
	ComplexValueType* dBuffer = this->m_Denominator->GetBufferPointer();
	size_t nPix = this->m_Denominator->GetLargestPossibleRegion().GetNumberOfPixels();
	for (size_t pix = 0; pix < nPix; pix++ ) {
		ComplexType cur = *(nBuffer+pix);
		*(nBuffer+pix) = cur * (*(dBuffer+pix));
	}
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::ComputeIterationEnergy() {
	VectorType* fBuffer = this->m_NextDeformationField->GetBufferPointer();
	size_t nPix = this->m_NextDeformationField->GetLargestPossibleRegion().GetNumberOfPixels();
	double totalNorm = 0;
	for (size_t pix = 0; pix < nPix; pix++ ) {
		totalNorm += (*(fBuffer+pix)).GetNorm();
	}

	this->m_CurrentLevelSetsValue = totalNorm/nPix;
}

}


#endif /* GRADIENTDESCENTLEVELSETSOPTIMIZER_HXX_ */
