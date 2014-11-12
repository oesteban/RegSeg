// --------------------------------------------------------------------------------------
// File:          OptimizerBase.hxx
// Date:          Feb 10, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
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
m_MaximumGradient(0.0),
m_MinimumConvergenceValue( 1e-3 ),
m_ConvergenceWindowSize( 10 ),
m_ConvergenceValue( 4.0 ),
m_Stop( false ),
m_StopCondition(MAXIMUM_NUMBER_OF_ITERATIONS),
m_CurrentIteration( 0 ),
m_NumberOfIterations( 250 ),
m_DescriptorRecomputationFreq(0),
m_ValueOscillations(0),
m_ValueOscillationsMax(1),
m_ValueOscillationsLast(0),
m_UseDescriptorRecomputation(false),
m_StepSize(1.0),
m_MaxSpeed(0.0),
m_MeanSpeed(0.0),
m_AvgSpeed(0.0),
m_AutoStepSize(true),
m_IsDiffeomorphic(true),
m_DiffeomorphismForced(false),
m_ForceDiffeomorphic(true),
m_UseLightWeightConvergenceChecking(true),
m_CurrentValue(itk::NumericTraits<MeasureType>::infinity()),
m_CurrentEnergy(itk::NumericTraits<MeasureType>::infinity()),
m_LastEnergy(itk::NumericTraits<MeasureType>::infinity()),
m_InitialValue(0.0)
{
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_GridSize.Fill( 0 );
}

template< typename TFunctional >
void OptimizerBase<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
	os << indent << "Learning rate:" << this->m_LearningRate << std::endl;
	os << indent << "Number of iterations: " << this->m_NumberOfIterations << std::endl;
	os << indent << "Current iteration: " << this->m_CurrentIteration << std::endl;
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

	if( this->m_UseDescriptorRecomputation ) {
		this->m_UseDescriptorRecomputation = (this->m_DescriptorRecomputationFreq > 0)
				&& (this->m_DescriptorRecomputationFreq < this->m_NumberOfIterations);
	}

	if( this->m_UseDescriptorRecomputation && this->m_ConvergenceWindowSize > this->m_DescriptorRecomputationFreq ) {
		itkWarningMacro( << "Convergence window size (" << this->m_ConvergenceWindowSize << ") is greater than descriptor recomputation frequency ("<<
				this->m_DescriptorRecomputationFreq << ")." << std::endl );
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
		if( this->m_CurrentIteration > this->m_DescriptorRecomputationFreq  && this->m_UseDescriptorRecomputation && ( (this->m_CurrentIteration -1) %this->m_DescriptorRecomputationFreq == 0 ) ) {
			this->m_Functional->UpdateDescriptors();
			// this->m_StepSize = this->m_StepSize * 1.0e-4;
			this->InvokeEvent( FunctionalModifiedEvent() );
			this->m_ConvergenceMonitoring->ClearEnergyValues();
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
			if (this->m_ConvergenceValue <= this->m_MinimumConvergenceValue) {
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

		if( this->m_MaximumGradient < 1e-5 ) {
			this->m_StopConditionDescription << "Maximum gradient update changed below the minimum threshold.";
			this->m_StopCondition = Self::STEP_TOO_SMALL;
			this->Stop();
			break;
		}


		if (this->m_InitialValue > 0.0){
			float inc = this->m_LastEnergy / this->m_CurrentEnergy;
			if (inc < 1.0) {
				this->m_ValueOscillations += 1;
				this->m_ValueOscillationsLast = this->m_CurrentIteration;
			} else if (inc > 1.0) {
				if ((this->m_CurrentIteration - this->m_ValueOscillationsLast) > this->m_ValueOscillationsMax) {
					float factor = 1.02 * inc * (1.0 - this->m_CurrentIteration/ this->m_NumberOfIterations);
					this->m_ValueOscillations = 0;
					this->m_StepSize *=  factor;
				}
			}

			if (this->m_ValueOscillations >= this->m_ValueOscillationsMax) {
				this->m_StepSize *= 0.5;
				this->m_ValueOscillations = 0;
				this->m_ConvergenceMonitoring->ClearEnergyValues();
			}
		} else {
			this->m_InitialValue = this->m_CurrentEnergy;
		}

		//std::cout << "StepSize=" << this->m_StepSize << std::endl;
		if( this->m_StepSize < 1e-8 ) {
			this->m_StopConditionDescription << "Parameters field changed below the minimum threshold.";
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
void OptimizerBase<TFunctional>
::SetStepSize (const InternalComputationValueType _arg) {
    if ( this->m_StepSize != _arg ){
    	this->m_AutoStepSize = false;
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
			("grid-size,S", bpo::value< size_t > (), "grid size")
			("update-descriptors,u", bpo::value< size_t > (), "frequency (iterations) to update descriptors of regions (0=no update)")
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
		this->SetDescriptorRecomputationFreq( updDesc );
	}

	bpo::variable_value v = this->m_Settings["convergence-energy"];
	this->m_UseLightWeightConvergenceChecking = ! v.as<bool>();
}

} // end namespace rstk




#endif /* OPTIMIZERBASE_HXX_ */
