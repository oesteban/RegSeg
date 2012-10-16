/* --------------------------------------------------------------------------------------
 * File:    EnergyOptimizer.h
 * Date:    15/10/2012
 * Author:  code@oscaresteban.es (Oscar Esteban)
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2012, Oscar Esteban - code@oscaresteban.es
 with Biomedical Image Technology, UPM (BIT-UPM) and
 Signal Processing Laboratory 5, EPFL (LTS5-EPFL).
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef EnergyOptimizer_h
#define EnergyOptimizer_h

#include <itkObject.h>
#include <itkObjectToObjectOptimizerBase.h>
#include "Energy.h"
#include <vector>

namespace rstk
{
/** \class EnergyOptimizer
 *  \bfief Abstract base class for variational optimizers.
 *
 *  \ingroup RSTK
 */
template< class TEnergy >
class EnergyOptimizer: public itk::ObjectToObjectOptimizerBase {
public:
	typedef EnergyOptimizer                            Self;
	typedef itk::ObjectToObjectOptimizerBase           Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( EnergyOptimizer, itk::ObjectToObjectOptimizerBase );

	/** Codes of stopping conditions. */
	typedef enum {
		MAXIMUM_NUMBER_OF_ITERATIONS,
		COSTFUNCTION_ERROR,
		UPDATE_PARAMETERS_ERROR,
		STEP_TOO_SMALL,
//		QUASI_NEWTON_STEP_ERROR,
		CONVERGENCE_CHECKER_PASSED,
		OTHER_ERROR
	} StopConditionType;

	itkNewMacro( Self );

	typedef TEnergy                                    EnergyType;
	typedef typename EnergyType::Pointer               EnergyPointer;
	typedef typename EnergyType::EnergyValueType       EnergyValueType;

	/** Derivative type */
//	typedef EnergyType::DerivativeType DerivativeType;
	typedef typename EnergyType::DisplacementFieldType DisplacementFieldType;
	typedef typename DisplacementFieldType::Pointer    DisplacementFieldPointer;

	/** Stop condition return string type */
	typedef std::string StopConditionReturnStringType;

	/** Stop condition internal string type */
	typedef std::ostringstream StopConditionDescriptionType;


	/** Measure type */
	typedef Superclass::MeasureType MeasureType;

	/** Internal computation type, for maintaining a desired precision */
	typedef Superclass::InternalComputationValueType InternalComputationValueType;

	/** Get the most recent gradient values. */
	itkGetConstReferenceMacro( Gradient, DerivativeType );

	/** Get stop condition enum */
	itkGetConstReferenceMacro(StopCondition, StopConditionType);

	/** Set the number of iterations. */
	itkSetMacro(NumberOfIterations, SizeValueType);

	/** Get the number of iterations. */
	itkGetConstReferenceMacro(NumberOfIterations, SizeValueType);

	/** Get the current iteration number. */
	itkGetConstMacro(CurrentIteration, SizeValueType);

	/** Resume optimization.
	 * This runs the optimization loop, and allows continuation
	 * of stopped optimization */
	virtual void ResumeOptimization() = 0;

	/** Stop optimization. The object is left in a state so the
	 * optimization can be resumed by calling ResumeOptimization. */
	virtual void StopOptimization(void);

	/** Get the reason for termination */
	virtual const StopConditionReturnStringType GetStopConditionDescription() const;


protected:
	EnergyOptimizer();
	virtual ~EnergyOptimizer() {}

	/* Common variables for optimization control and reporting */
	bool                          m_Stop;
	StopConditionType             m_StopCondition;
	StopConditionDescriptionType  m_StopConditionDescription;
	SizeValueType                 m_NumberOfIterations;
	SizeValueType                 m_CurrentIteration;

	/** Current gradient */
	EnergyPointer m_Energy;
	DerivativeType     m_Gradient;
	virtual void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
	EnergyOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "EnergyOptimizer.txx"
#endif


#endif /* EnergyOptimizer_h */
