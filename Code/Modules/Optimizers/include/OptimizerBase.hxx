// This file is part of RegSeg
//
// Copyright 2014-2017, Oscar Esteban <code@oscaresteban.es>
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#ifndef OPTIMIZERBASE_HXX_
#define OPTIMIZERBASE_HXX_


#include "OptimizerBase.h"

#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <itkComplexToRealImageFilter.h>
#include <itkImageAlgorithm.h>

using namespace std;

#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

namespace rstk {

/**
 * Default constructor
 */
template< typename TFunctional >
OptimizerBase<TFunctional>::OptimizerBase():
m_LearningRate( 1.0 ),
m_Momentum( 0.0 ),
m_LastMaximumGradient(0.0),
m_MaximumGradient(itk::NumericTraits<InternalComputationValueType>::infinity()),
m_MinimumConvergenceValue( 1e-5 ),
m_ConvergenceWindowSize( 10 ),
m_ConvergenceValue( itk::NumericTraits<InternalComputationValueType>::infinity() ),
m_Stop( false ),
m_StopCondition(MAXIMUM_NUMBER_OF_ITERATIONS),
m_DescriptorRecompPeriod(0),
m_NextRecompIteration(1),
m_ValueOscillations(0),
m_ValueOscillationsMax(1),
m_ValueOscillationsLast(0),
m_UseDescriptorRecomputation(false),
m_StepSize(1.0),
m_MaxSpeed(0.0),
m_MeanSpeed(0.0),
m_AvgSpeed(0.0),
m_AutoStepSize(false),
m_IsDiffeomorphic(true),
m_DiffeomorphismForced(false),
m_ForceDiffeomorphic(true),
m_UseLightWeightConvergenceChecking(true),
m_UseAdaptativeDescriptors(false),
m_CurrentValue(itk::NumericTraits<MeasureType>::infinity()),
m_CurrentEnergy(itk::NumericTraits<MeasureType>::infinity()),
m_CurrentNorm(0.0),
m_LastEnergy(itk::NumericTraits<MeasureType>::infinity())
{
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_GridSize.Fill( 0 );
	this->m_GridSpacing.Fill(0.0);
}

template< typename TFunctional >
void OptimizerBase<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
	os << indent << "Learning rate:" << this->m_LearningRate << std::endl;
	os << indent << "Stop condition:" << this->m_StopCondition << std::endl;
	os << indent << "Stop condition description: " << this->m_StopConditionDescription.str() << std::endl;
}


template< typename TFunctional >
void OptimizerBase<TFunctional>::Start() {
	itkDebugMacro("OptimizerBase::Start()");

	if ( this->m_Settings.size() > 0 ) {
		this->ParseSettings();
	}

	/* Settings validation */
	if ( this->m_Functional.IsNull() ) {
		itkExceptionMacro("Energy functional must be set");
	}

	this->m_Functional->Initialize();

	/* Check & initialize parameter fields */
	this->InitializeParameters();
	this->InitializeAuxiliarParameters();

	if (this->m_UseAdaptativeDescriptors ) {
		this->m_UseDescriptorRecomputation = true;
		this->m_DescriptorRecompPeriod = 1;
	}

	if( this->m_UseDescriptorRecomputation ) {
		this->m_UseDescriptorRecomputation = (this->m_DescriptorRecompPeriod > 0)
				&& (this->m_DescriptorRecompPeriod < this->m_NumberOfIterations);
	}

	/* Initialize convergence checker */
	this->m_ConvergenceMonitoring = ConvergenceMonitoringType::New();
	this->m_ConvergenceMonitoring->SetWindowSize( this->m_ConvergenceWindowSize );

//	if( this->m_ReturnBestParametersAndValue )	{
//		this->m_BestParameters = this->GetCurrentPosition( );
//		this->m_CurrentBestValue = NumericTraits< MeasureType >::max();
//	}

	this->InvokeEvent( itk::StartEvent() );
	this->m_CurrentIteration++;
	this->Resume();
}

template< typename TFunctional >
const typename OptimizerBase<TFunctional>::StopConditionReturnStringType
OptimizerBase<TFunctional>::
GetStopConditionDescription() const {
  return this->m_StopConditionDescription.str();
}

template< typename TFunctional >
void OptimizerBase<TFunctional>::Stop() {
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
void OptimizerBase<TFunctional>::Resume() {
	this->m_StopConditionDescription.str("");
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_Stop = false;

	while( ! this->m_Stop )	{
		if( this->DoDescriptorsUpdate()) {
			this->m_Functional->UpdateDescriptors();
			this->InvokeEvent( FunctionalModifiedEvent() );
		}


		/* Compute functional value/derivative. */
		try	{
			this->ComputeDerivative();
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

		/* Update the deformation field */
		this->PostIteration();


		/* TODO Store best value and position */
		//if ( this->m_ReturnBestParametersAndValue && this->m_CurrentValue < this->m_CurrentBestValue )
		//{
		//	this->m_CurrentBestValue = this->m_CurrentValue;
		//	this->m_BestParameters = this->GetCurrentPosition( );
		//}


		/*
		 * Check the convergence by WindowConvergenceMonitoringFunction.
		 */
		this->m_ConvergenceMonitoring->AddEnergyValue( this->m_CurrentValue );

		try {
			this->m_ConvergenceValue = this->m_ConvergenceMonitoring->GetConvergenceValue();
			if (fabs(this->m_ConvergenceValue) <= this->m_MinimumConvergenceValue) {
				this->m_StopConditionDescription << "Convergence checker passed at iteration " << this->m_CurrentIteration << ".";
				this->m_StopCondition = Self::CONVERGENCE_CHECKER_PASSED;
				this->Stop();
				break;
			}
		}
		catch(std::exception & e) {
			std::cerr << "GetConvergenceValue() failed with exception: " << e.what() << std::endl;
		}

		/* Update and check iteration count */
		if ( this->m_CurrentIteration >= this->m_NumberOfIterations ) {
			this->m_StopConditionDescription << "Maximum number of iterations (" << this->m_NumberOfIterations << ") exceeded.";
			this->m_StopCondition = MAXIMUM_NUMBER_OF_ITERATIONS;
			this->Stop();
			break;
		}

		if( (this->m_MaximumGradient * this->m_StepSize) < 1e-5 ) {
			this->m_StopConditionDescription << "Maximum gradient update changed below the minimum threshold.";
			this->m_StopCondition = Self::STEP_TOO_SMALL;
			this->Stop();
			break;
		}


		if (this->m_CurrentIteration > 0){
			float inc = 1.0;
			if (this->m_ConvergenceValue != itk::NumericTraits<InternalComputationValueType>::infinity() ) {
				inc+= this->m_ConvergenceValue;
			} else {
				inc = this->m_LastEnergy / this->m_CurrentEnergy;
			}
			if (inc < 1.0) {
				this->m_ValueOscillations += 1;
				this->m_ValueOscillationsLast = this->m_CurrentIteration;
			} else if (inc >= 1.0) {
				if ((this->m_CurrentIteration - this->m_ValueOscillationsLast) > this->m_ValueOscillationsMax) {
					float factor = inc * (1.0 - ((this->m_CurrentIteration * this->m_CurrentIteration) / this->m_NumberOfIterations));
					if (factor < 1.0) {
						factor = 1.0;
					}
					this->m_ValueOscillations = 0;
					this->m_StepSize *=  factor;
				}
			}

			if (this->m_ValueOscillations >= this->m_ValueOscillationsMax) {
				this->m_StepSize *= 0.5;
				this->m_ValueOscillations = 0;
			}
		}

		if( this->m_StepSize < 1e-8 ) {
			this->m_StopConditionDescription << "step size fell below the minimum.";
			this->m_StopCondition = Self::STEP_TOO_SMALL;
			this->Stop();
			break;
		}
		this->m_LastEnergy = this->m_CurrentEnergy;
		this->m_LastMaximumGradient = this->m_MaximumGradient;

		this->InvokeEvent( itk::IterationEvent() );
		this->m_CurrentIteration++;
	} //while (!m_Stop)
}

template< typename TFunctional >
bool OptimizerBase<TFunctional>
::DoDescriptorsUpdate() {
	if (!this->m_UseDescriptorRecomputation ) return false;

	size_t lastIt = this->m_CurrentIteration -1;
	if ( !this->m_UseAdaptativeDescriptors ) return lastIt > 0 && (lastIt % this->m_DescriptorRecompPeriod) == 0;

	bool val = false;
	if( this->m_CurrentIteration == m_NextRecompIteration ) {
		this->m_NextRecompIteration+= int(ceil( 1.0/exp(-(this->m_CurrentIteration * 1.5 / this->m_ConvergenceWindowSize ))));
		val = true;
	}
	this->m_UseDescriptorRecomputation = this->m_CurrentIteration < this->m_ConvergenceWindowSize;
	return val;
	// this->m_ConvergenceMonitoring->ClearEnergyValues();
}


template< typename TFunctional >
void OptimizerBase<TFunctional>
::SetStepSize (const InternalComputationValueType _arg) {
    if ( this->m_StepSize != _arg ){
    	this->m_StepSize = _arg;
    	this->Modified();
    }
}

template< typename TFunctional >
void OptimizerBase<TFunctional>
::AddOptions( SettingsDesc& opts ) {
	opts.add_options()
			("alpha,a", bpo::value< std::vector<float> >()->multitoken(), "alpha value in regularization")
			("beta,b", bpo::value< std::vector<float> >()->multitoken(), "beta value in regularization")
			("step-size,s", bpo::value< double > (), "step-size value in optimization")
			("gradient-scales,g", bpo::value< std::vector<float> >()->multitoken(), "alpha value in regularization")
			("learning-rate,r", bpo::value< float > (), "learning rate to update step size")
			("iterations,i", bpo::value< size_t > (), "number of iterations")
			("convergence-window,w", bpo::value< size_t > (), "number of iterations of convergence window")
			("convergence-thresh,t", bpo::value< double > (), "convergence value")
			("grid-size", bpo::value< std::vector<size_t> >()->multitoken(), "size of control points grid")
			("grid-spacing", bpo::value< std::vector<float> >()->multitoken(), "spacing between control points ")
			("update-descriptors,u", bpo::value< size_t > (), "frequency (iterations) to update descriptors of regions (0=no update)")
			("adaptative-descriptors", bpo::bool_switch(), "recomputes descriptors more often at the beginning of the process")
			("step-auto", bpo::bool_switch(), "guess appropriate step size depending on first iteration")
			("convergence-energy", bpo::bool_switch(), "disables lazy convergence tracking: instead of fast computation of the mean norm of "
					"the displacement field, it computes the full energy functional");
}

template< typename TFunctional >
void OptimizerBase<TFunctional>
::ParseSettings() {
	bpo::notify( this->m_Settings );

	if( this->m_Settings.count( "step-size" ) ){
		bpo::variable_value v = this->m_Settings["step-size"];
		this->SetStepSize( v.as< double >() );
	}

	if( this->m_Settings.count( "gradient-scales" ) ){
		bpo::variable_value v = this->m_Settings["gradient-scales"];
		std::vector<float> s = v.as< std::vector<float> > ();
		ScalesType scales(Dimension);
		if (scales.size() == 1) {
			scales.Fill(s[0]);
		} else if (scales.size() == Dimension) {
			for( size_t i = 0; i < Dimension; i++)
				scales[i] = s[i];
		}
		this->SetScales(scales);
	}

	if( this->m_Settings.count( "learning-rate" ) ){
		bpo::variable_value v = this->m_Settings["learning-rate"];
		this->SetLearningRate( v.as< float >() );
	}

	if( this->m_Settings.count( "iterations" ) ){
		bpo::variable_value v = this->m_Settings["iterations"];
		this->SetNumberOfIterations( v.as< size_t >() );

		if( ! this->m_Settings.count("convergence-window" ) ) {
			this->m_ConvergenceWindowSize = floor( 0.20 * this->m_NumberOfIterations ) + 1;
		}
	}

	if( this->m_Settings.count( "convergence-window" ) ){
			bpo::variable_value v = this->m_Settings["convergence-window"];
			this->m_ConvergenceWindowSize = v.as< size_t >();
	}

	if( this->m_Settings.count( "convergence-thresh" ) ){
			bpo::variable_value v = this->m_Settings["convergence-thresh"];
			this->m_MinimumConvergenceValue = v.as< double >();
	}

	if (this->m_Settings.count("update-descriptors")) {
		bpo::variable_value v = this->m_Settings["update-descriptors"];
		size_t updDesc =  v.as<size_t>();
		this->SetUseDescriptorRecomputation(true);
		this->SetDescriptorRecompPeriod( updDesc );
	}

	bpo::variable_value v = this->m_Settings["convergence-energy"];
	this->m_UseLightWeightConvergenceChecking = ! v.as<bool>();

	bpo::variable_value a = this->m_Settings["adaptative-descriptors"];
	this->m_UseAdaptativeDescriptors = a.as<bool>();

	bpo::variable_value sa = this->m_Settings["step-auto"];
	this->m_AutoStepSize = sa.as<bool>();
}

} // end namespace rstk




#endif /* OPTIMIZERBASE_HXX_ */
