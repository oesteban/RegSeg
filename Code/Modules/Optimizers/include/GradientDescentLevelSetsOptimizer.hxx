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
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <cmath>
#include <itkComplexToRealImageFilter.h>
#include "DisplacementFieldFileWriter.h"
#include <itkMeshFileWriter.h>
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;

namespace rstk {

/**
 * Default constructor
 */
template< typename TLevelSetsFunction >
GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::GradientDescentLevelSetsOptimizer() {
	this->m_LearningRate = itk::NumericTraits<InternalComputationValueType>::One;
	this->m_MaximumStepSizeInPhysicalUnits = itk::NumericTraits<InternalComputationValueType>::Zero;
	this->m_MinimumConvergenceValue = 1e-8;
	this->m_ConvergenceWindowSize = 30;
	this->m_StepSize = 0.5;
	this->m_A.SetIdentity();
	this->m_A(0,0) = 10e3;
	this->m_A(1,1) = 10e3;
	this->m_A(2,2) = 10e10;
	this->m_B.SetIdentity();
	this->m_B*= 10e4;
	this->m_B(2,2) = 10e10;
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
		PointType orig = this->m_LevelSetsFunction->GetReferenceImage()->GetOrigin();
		PointType end;

		typename DeformationFieldType::IndexType end_idx;
		typename DeformationFieldType::SizeType refSize = this->m_LevelSetsFunction->GetReferenceImage()->GetLargestPossibleRegion().GetSize();
		for( size_t dim = 0; dim < DeformationFieldType::ImageDimension; dim++ ) end_idx[dim] = refSize[dim];
		this->m_LevelSetsFunction->GetReferenceImage()->TransformIndexToPhysicalPoint( end_idx, end );
		this->InitializeDeformationField( orig, end, this->m_LevelSetsFunction->GetReferenceImage()->GetDirection() );
	}

	/* Initialize next deformationfield */
	this->m_NextDeformationField = DeformationFieldType::New();
	this->m_NextDeformationField->SetRegions( this->m_DeformationField->GetLargestPossibleRegion() );
	this->m_NextDeformationField->SetSpacing( this->m_DeformationField->GetSpacing() );
	this->m_NextDeformationField->SetDirection( this->m_DeformationField->GetDirection() );
	this->m_NextDeformationField->SetOrigin( this->m_DeformationField->GetOrigin() );
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

	this->m_LevelSetsFunction->CopyInformation( this->m_DeformationField );
	this->m_LevelSetsFunction->Initialize();

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
			this->m_LevelSetsFunction->ComputeGradient();
#ifndef NDEBUG
			typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;
			typename Writer::Pointer p = Writer::New();
			std::stringstream ss2;
			ss2 << "gradient_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
			p->SetFileName( ss2.str().c_str() );
			p->SetInput( this->m_LevelSetsFunction->GetGradientMap() );
			p->Update();
#endif
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

		this->ComputeIterationChange();

#ifndef NDEBUG
			typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;
			typename Writer::Pointer p = Writer::New();
			std::stringstream ss2;
			ss2 << "field_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
			p->SetFileName( ss2.str().c_str() );
			p->SetInput( this->m_NextDeformationField );
			p->Update();
#endif

		/* Update the level sets contour and deformation field */
		double updateNorm = this->m_LevelSetsFunction->UpdateContour( this->m_NextDeformationField );

		const VectorType* next = this->m_NextDeformationField->GetBufferPointer();
		VectorType* curr = this->m_DeformationField->GetBufferPointer();
		size_t nPix = this->m_DeformationField->GetLargestPossibleRegion().GetNumberOfPixels();

		for( size_t pix = 0; pix<nPix; pix++){
			*(curr+pix) = *(next+pix);
		}


		/* TODO Store best value and position */
		//if ( this->m_ReturnBestParametersAndValue && this->m_CurrentLevelSetsValue < this->m_CurrentBestValue )
		//{
		//	this->m_CurrentBestValue = this->m_CurrentLevelSetsValue;
		//	this->m_BestParameters = this->GetCurrentPosition( );
		//}

#ifndef NDEBUG
			typedef itk::MeshFileWriter< ContourType >     MeshWriterType;
			typename MeshWriterType::Pointer w = MeshWriterType::New();

			for ( size_t c = 0; c < this->m_LevelSetsFunction->GetCurrentContourPosition().size(); c++ ) {
				w->SetInput( this->m_LevelSetsFunction->GetCurrentContourPosition()[c] );
				std::stringstream ss;
				ss << "contour_" << std::setfill('0') << std::setw(2) << c <<"_it_" << std::setw(3) << this->m_CurrentIteration << ".vtk";
				w->SetFileName( ss.str() );
				w->Update();
			}

			typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;
			typename Writer::Pointer p = Writer::New();
			std::stringstream ss2;
			ss2 << "field_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
			p->SetFileName( ss2.str().c_str() );
			p->SetInput( this->m_DeformationField );
			p->Update();
#endif

		std::cout << "[" << this->m_CurrentIteration << "] " << this->m_CurrentLevelSetsValue << " | " << updateNorm << " | " << this->m_LevelSetsFunction->GetValue() << std::endl;

		/*
		 * Check the convergence by WindowConvergenceMonitoringFunction.
		 */

		/*
		this->m_ConvergenceMonitoring->AddEnergyValue( this->m_CurrentLevelSetsValue );
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


		if( this->m_CurrentLevelSetsValue < this->m_ConvergenceValue ) {
			this->m_StopConditionDescription << "Deformation field changed below the minimum threshold.";
			this->m_StopCondition = Superclass::STEP_TOO_SMALL;
			this->Stop();
			break;
		}
		*/

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
	double min_val = 1.0e-5;

	DeformationFieldConstPointer shapeGradient = this->m_LevelSetsFunction->GetGradientMap();

#ifndef NDEBUG
	typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;
	typename Writer::Pointer p = Writer::New();
	std::stringstream ss2;
	ss2 << "gradient2_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
	p->SetFileName( ss2.str().c_str() );
	p->SetInput( shapeGradient );
	p->Update();

	typename Writer::Pointer p2 = Writer::New();
	ss2.str("");
	ss2 << "field2_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
	p->SetFileName( ss2.str().c_str() );
	p->SetInput( shapeGradient );
	p->Update();
#endif

	// Get deformation fields pointers
	const VectorType* dFieldBuffer = this->m_DeformationField->GetBufferPointer();
	const VectorType* sFieldBuffer = shapeGradient->GetBufferPointer();

	// Create container for deformation field components
	DeformationComponentPointer fieldComponent = DeformationComponentType::New();
	fieldComponent->SetDirection( this->m_DeformationField->GetDirection() );
	fieldComponent->SetRegions( this->m_DeformationField->GetLargestPossibleRegion() );
	fieldComponent->Allocate();
	PointValueType* dCompBuffer = fieldComponent->GetBufferPointer();

	ComplexFieldPointer fftNum;

	size_t nPix = fieldComponent->GetLargestPossibleRegion().GetNumberOfPixels();
	for( size_t d = 0; d < Dimension; d++ ) {
		// Get component from vector image
		size_t comp, idx2;
		PointValueType u,grad;
		for(size_t pix = 0; pix< nPix; pix++) {
			u = (*(dFieldBuffer+pix))[d];
			grad = (*(sFieldBuffer+pix))[d];
			*(dCompBuffer+pix) = u / this->m_StepSize + grad; // ATTENTION: the + symbol depends on the definition of the normal
		}                                                                                                                 // OUTWARDS-> + ; INWARDS-> -.

		FFTPointer fftFilter = FFTType::New();
		fftFilter->SetInput( fieldComponent );
		fftFilter->Update();  // This is required for computing the denominator for first time
		FTDomainPointer fftComponent =  fftFilter->GetOutput();

		if ( fftNum.IsNull() ) {
			fftNum = ComplexFieldType::New();
			fftNum->SetRegions( fftComponent->GetLargestPossibleRegion() );
			//fftNum->CopyInformation( fftComponent );
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


	VectorType* nextFieldBuffer = this->m_NextDeformationField->GetBufferPointer();
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
			this->m_StopCondition = Superclass::UPDATE_PARAMETERS_ERROR;
			this->m_StopConditionDescription << "UpdateDeformationField error";
			this->Stop();
			// Pass exception to caller
			throw err;
		}

		// Set component on this->m_NextDeformationField
		const PointValueType* resultBuffer = ifft->GetOutput()->GetBufferPointer();
		PointValueType val;
		for(size_t pix = 0; pix< nPix; pix++) {
			val = *(resultBuffer+pix);
			(*(nextFieldBuffer+pix))[d] = ( fabs(val) > min_val )?val:0.0;
		}
	}

	this->InvokeEvent( itk::IterationEvent() );
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>
::ApplyRegularizationTerm( typename GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::ComplexFieldType* reference ){
	if( this->m_Denominator.IsNull() ) { // If executed for first time, cache the map
		this->InitializeDenominator( reference );
	}

	// Fill reference buffer with elements updated with the filter
	ComplexFieldValue* nBuffer = reference->GetBufferPointer();
	MatrixType* dBuffer = this->m_Denominator->GetBufferPointer();
	size_t nPix = this->m_Denominator->GetLargestPossibleRegion().GetNumberOfPixels();
	ComplexFieldValue curval;
	ComplexValuesVector realval;
	ComplexValuesVector imagval;

	for (size_t pix = 0; pix < nPix; pix++ ) {
		curval = *(nBuffer+pix);
		MatrixType tensor = *(dBuffer+pix);

		for( size_t i = 0; i<Dimension; i++ ) {
			realval[i] = curval[i].real();
			imagval[i] = curval[i].imag();
		}

		realval = tensor * realval;
		imagval = tensor * imagval;

		for( size_t i = 0; i<Dimension; i++ ) {
			curval[i] = ComplexType(realval[i],imagval[i]);
		}


		*(nBuffer+pix) = curval;
	}
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>
::InitializeDenominator( typename GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::ComplexFieldType* reference ){
	double pi2 = 2* vnl_math::pi;
	MatrixType I, ddor, zero;
	I.SetIdentity();
	zero.Fill(0.0);
	MatrixType C = ( I * (1.0/this->m_StepSize) + m_A*I ) * pi2;

	this->m_Denominator = TensorFieldType::New();
	this->m_Denominator->SetSpacing(   reference->GetSpacing() );
	this->m_Denominator->SetDirection( reference->GetDirection() );
	this->m_Denominator->SetRegions( reference->GetLargestPossibleRegion().GetSize() );
	this->m_Denominator->Allocate();
	this->m_Denominator->FillBuffer( zero );

	typename TensorFieldType::SizeType size = this->m_Denominator->GetLargestPossibleRegion().GetSize();
	MatrixType* buffer = this->m_Denominator->GetBufferPointer();
	size_t nPix = this->m_Denominator->GetLargestPossibleRegion().GetNumberOfPixels();

	// Fill Buffer with denominator data
	double lag_el; // accumulates the FT{lagrange operator}.
	typename TensorFieldType::IndexType idx;
	for (size_t pix = 0; pix < nPix; pix++ ) {
		lag_el = 0.0;
		idx = this->m_Denominator->ComputeIndex( pix );
		for(size_t d = 0; d < Dimension; d++ )
			lag_el+= 2.0 * cos( (pi2*idx[d])/size[d]) - 2.0;

		ddor = C - m_B * lag_el;
		*(buffer+pix) = ddor.GetInverse();
	}
}

template< typename TLevelSetsFunction >
void GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::ComputeIterationChange() {
	VectorType* fnextBuffer = this->m_NextDeformationField->GetBufferPointer();
	VectorType* fBuffer = this->m_DeformationField->GetBufferPointer();
	size_t nPix = this->m_NextDeformationField->GetLargestPossibleRegion().GetNumberOfPixels();
	double totalNorm = 0;
	VectorType t0,t1;
	for (size_t pix = 0; pix < nPix; pix++ ) {
		t0 = *(fBuffer+pix);
		t1 = *(fnextBuffer+pix);
		totalNorm += ( t1 - t0 ).GetNorm();
	}

	this->m_CurrentLevelSetsValue = totalNorm/nPix;
}


template< typename TLevelSetsFunction >
void
GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::InitializeDeformationField(
	 const typename GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::PointType orig,
	 const typename GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::PointType end,
	 const typename GradientDescentLevelSetsOptimizer<TLevelSetsFunction>::DeformationFieldDirectionType dir)
{
	this->m_DeformationField = DeformationFieldType::New();

	typename DeformationFieldType::SpacingType spacing;
	for( size_t dim = 0; dim < DeformationFieldType::ImageDimension; dim++ ) {
		spacing[dim] = fabs( (end[dim]-orig[dim])/(1.0*this->m_GridSize[dim]));
	}

	VectorType zerov;
	zerov.Fill(0.0);

	this->m_DeformationField->SetRegions( this->m_GridSize );
	this->m_DeformationField->SetSpacing( spacing );
	this->m_DeformationField->SetDirection( dir );
	this->m_DeformationField->SetOrigin( orig );
	this->m_DeformationField->Allocate();
	this->m_DeformationField->FillBuffer( zerov );
}

}
#endif /* GRADIENTDESCENTLEVELSETSOPTIMIZER_HXX_ */
