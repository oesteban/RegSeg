/* --------------------------------------------------------------------------------------
 * File:    GradientDescentEnergyOptimizer.h
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

#ifndef GradientDescentEnergyOptimizer_h
#define GradientDescentEnergyOptimizer_h

#include <itkObject.h>
#include "EnergyOptimizer.h"
#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>

namespace rstk
{

template< class TEnergy >
class GradientDescentEnergyOptimizer: public EnergyOptimizer< TEnergy > {
public:
	/** Standard class typedefs and macros */
	typedef GradientDescentEnergyOptimizer             Self;
	typedef EnergyOptimizer< TEnergy >                 Superclass;
	typedef SmartPointer<Self>                         Pointer;
	typedef SmartPointer< const Self >                 ConstPointer;
	itkTypeMacro( GradientDescentEnergyOptimizer, EnergyOptimizer ); // Run-time type information (and related methods)
	itkNewMacro( Self );                                             // New macro for creation of through a Smart Pointer

	/** Energy-related typedefs */
	typedef typename Superclass::EnergyType            EnergyType;
	typedef typename EnergyType::Pointer               EnergyPointer;
	typedef typename EnergyType::EnergyValueType       EnergyValueType;
	typedef typename EnergyType::DisplacementFieldType DisplacementFieldType;
	typedef typename DisplacementFieldType::Pointer    DisplacementFieldPointer;

	/** Derivative type */
	typedef Superclass::DerivativeType      DerivativeType;

	/** Metric type over which this class is templated */
	typedef Superclass::MeasureType                  MeasureType;
	typedef Superclass::InternalComputationValueType InternalComputationValueType;

	/** Type for the convergence checker */
	typedef itk::Function::WindowConvergenceMonitoringFunction<double>	ConvergenceMonitoringType;

	itkSetMacro(LearningRate, InternalComputationValueType);               // Set the learning rate
	itkGetConstReferenceMacro(LearningRate, InternalComputationValueType); // Get the learning rate

	/** Minimum convergence value for convergence checking.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the energy profile. When the convergence value reaches
	 *  a small value, it would be treated as converged.
	 *
	 *  The default m_MinimumConvergenceValue is set to 1e-8 to pass all
	 *  tests. It is suggested to use 1e-6 for less stringent convergence
	 *  checking.
	 */
	itkSetMacro(MinimumConvergenceValue, InternalComputationValueType);

	/** Window size for the convergence checker.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the energy (metric value) profile.
	 *
	 *  The default m_ConvergenceWindowSize is set to 50 to pass all
	 *  tests. It is suggested to use 10 for less stringent convergence
	 *  checking.
	 */
	itkSetMacro(ConvergenceWindowSize, SizeValueType);

	/** Get current convergence value */
	itkGetConstReferenceMacro( ConvergenceValue, InternalComputationValueType );

	/** Flag. Set to have the optimizer track and return the best
	 *  best metric value and corresponding best parameters that were
	 *  calculated during the optimization. This captures the best
	 *  solution when the optimizer oversteps or osciallates near the end
	 *  of an optimization.
	 *  Results are stored in m_CurrentMetricValue and in the assigned metric's
	 *  parameters, retrievable via optimizer->GetCurrentPosition().
	 *  This option requires additional memory to store the best
	 *  parameters, which can be large when working with high-dimensional
	 *  transforms such as DisplacementFieldTransform.
	 */
	itkSetMacro(ReturnBestParametersAndValue, bool);
	itkGetConstReferenceMacro(ReturnBestParametersAndValue, bool);
	itkBooleanMacro(ReturnBestParametersAndValue);

	/** Start and run the optimization */
	virtual void StartOptimization();

	virtual void StopOptimization(void);

	virtual void ResumeOptimization();

protected:
	GradientDescentEnergyOptimizer();
	virtual ~GradientDescentEnergyOptimizer() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const {
		Superclass::PrintSelf( os, indent );
	}

private:
	EnergyPointer m_Energy;

	GradientDescentEnergyOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "GradientDescentEnergyOptimizer.txx"
#endif


#endif /* GradientDescentEnergyOptimizer_h */
