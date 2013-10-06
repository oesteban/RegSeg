// --------------------------------------------------------------------------------------
// File:          SpectralOptimizer.h
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

#ifndef SPECTRALOPTIMIZER_H_
#define SPECTRALOPTIMIZER_H_


#include <itkObjectToObjectOptimizerBase.h>

#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkRealToHalfHermitianForwardFFTImageFilter.h>
#include <itkHalfHermitianToRealInverseFFTImageFilter.h>

#include <itkImageIteratorWithIndex.h>
#include <itkImageAlgorithm.h>

using namespace itk;

namespace rstk
{
/**
 * \class SpectralOptimizer
 *  \brief Gradient descent optimizer.
 *
 * GradientDescentOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current deformation field is updated according:
 * \f[
 *        u^{t+1} = \mathcal{FT^{-1}}
 * \f]
 */

template< typename TFunctional >
class SpectralOptimizer: public ObjectToObjectOptimizerBase {
public:
	/** Standard class typedefs and macros */
	typedef SpectralOptimizer           Self;
	typedef ObjectToObjectOptimizerBase                Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( SpectralOptimizer, ObjectToObjectOptimizerBase ); // Run-time type information (and related methods)

	/** Metric type over which this class is templated */
	typedef TFunctional                                             FunctionalType;
	itkStaticConstMacro( Dimension, unsigned int, FunctionalType::Dimension );

	/** Codes of stopping conditions. */
	typedef enum {
		MAXIMUM_NUMBER_OF_ITERATIONS,
		COSTFUNCTION_ERROR,
		UPDATE_PARAMETERS_ERROR,
		STEP_TOO_SMALL,
		QUASI_NEWTON_STEP_ERROR,
		CONVERGENCE_CHECKER_PASSED,
		OTHER_ERROR
	} StopConditionType;

	/** Stop condition return string type */
	typedef std::string                                             StopConditionReturnStringType;

	/** Stop condition internal string type */
	typedef std::ostringstream                                      StopConditionDescriptionType;

	/** Internal computation value type */
	typedef double                                                  InternalComputationValueType;

	/** Functional definitions */
	typedef typename FunctionalType::Pointer                        FunctionalPointer;
	typedef typename FunctionalType::FieldType                      ParametersType;
	typedef typename FunctionalType::FieldType                      DerivativeType;
	typedef typename FunctionalType::MeasureType                    MeasureType;

	typedef typename FunctionalType::PointType                      PointType;
	typedef typename FunctionalType::VectorType                     VectorType;
	typedef typename FunctionalType::PointValueType                 PointValueType;

	typedef typename ParametersType::Pointer                        ParametersPointer;
	typedef typename ParametersType::ConstPointer                   ParametersConstPointer;
	typedef typename ParametersType::PointType                      ParametersPointType;
	typedef typename ParametersType::DirectionType                  ParametersDirectionType;
	typedef typename ParametersType::SizeType                       GridSizeType;
	typedef typename itk::Image<PointValueType, Dimension >         ParametersComponentType;
	typedef typename ParametersComponentType::Pointer               ParametersComponentPointer;

	typedef itk::RealToHalfHermitianForwardFFTImageFilter
			                          <ParametersComponentType>     FFTType;
	typedef typename FFTType::Pointer                               FFTPointer;
	typedef typename FFTType::OutputImageType                       FTDomainType;
	typedef typename FTDomainType::Pointer                          FTDomainPointer;
	typedef typename FTDomainType::PixelType                        ComplexType;
	typedef typename ComplexType::value_type                        ComplexValueType;

	typedef itk::HalfHermitianToRealInverseFFTImageFilter
			        <FTDomainType, ParametersComponentType>         IFFTType;
	typedef typename IFFTType::Pointer                              IFFTPointer;


	typedef itk::Image< ComplexValueType, Dimension >               RealPartType;
	typedef itk::Vector< ComplexValueType, Dimension >              ComplexValuesVector;

	typedef itk::Vector< ComplexType, Dimension >                   ComplexFieldValue;
	typedef itk::Image< ComplexFieldValue, Dimension >              ComplexFieldType;
	typedef typename ComplexFieldType::Pointer                      ComplexFieldPointer;

	typedef itk::Matrix< ComplexValueType, Dimension, Dimension >   MatrixType;
	typedef itk::Image< MatrixType, Dimension >                     TensorFieldType;
	typedef typename TensorFieldType::Pointer                       TensorFieldPointer;

	typedef size_t SizeValueType;

	/** Type for the convergence checker */
	typedef itk::Function::WindowConvergenceMonitoringFunction<MeasureType>	         ConvergenceMonitoringType;

	/** Accessors for Functional */
	itkGetObjectMacro( Functional, FunctionalType );
	itkSetObjectMacro( Functional, FunctionalType );

	itkSetMacro(LearningRate, InternalComputationValueType);               // Set the learning rate
	itkGetConstReferenceMacro(LearningRate, InternalComputationValueType); // Get the learning rate

	itkSetObjectMacro(Parameters, ParametersType);
	itkGetConstObjectMacro(Parameters, ParametersType);

	/** Minimum convergence value for convergence checking.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the Functional profile. When the convergence value reaches
	 *  a small value, it would be treated as converged.
	 *
	 *  The default m_MinimumConvergenceValue is set to 1e-8 to pass all
	 *  tests. It is suggested to use 1e-6 for less stringent convergence
	 *  checking.
	 */
	itkSetMacro(MinimumConvergenceValue, InternalComputationValueType);

	/** Window size for the convergence checker.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the Functional (metric value) profile.
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
//	itkSetMacro(ReturnBestParametersAndValue, bool);
//	itkGetConstReferenceMacro(ReturnBestParametersAndValue, bool);
//	itkBooleanMacro(ReturnBestParametersAndValue);

	itkSetMacro( A, MatrixType );
	itkGetConstMacro( A, MatrixType );

	itkSetMacro( B, MatrixType );
	itkGetConstMacro( B, MatrixType );

	itkSetMacro( StepSize, InternalComputationValueType );
	itkGetConstMacro( StepSize, InternalComputationValueType );

	void SetAlpha( InternalComputationValueType v ) {
		this->m_A.SetIdentity();
		this->m_A = this->m_A * 2.0 * v;
	}

	void SetBeta( InternalComputationValueType v ) {
		this->m_B.SetIdentity();
		this->m_B = this->m_B * 2.0 * v;
	}

	/** Get stop condition enum */
	itkGetConstReferenceMacro(StopCondition, StopConditionType);

	/** Set the number of iterations. */
	itkSetMacro(NumberOfIterations, SizeValueType);

	itkSetMacro( GridSize, GridSizeType );

	void SetGridSize( double val ) { this->m_GridSize.Fill(val); }

	/** Get the number of iterations. */
	itkGetConstReferenceMacro(NumberOfIterations, SizeValueType);

	/** Get the current iteration number. */
	itkGetConstMacro(CurrentIteration, SizeValueType);

	/** Start and run the optimization */
	void Start();

	void Stop(void);

	/** Get the reason for termination */
	const StopConditionReturnStringType GetStopConditionDescription() const;

	void Resume();


	MeasureType ComputeIterationChange();
	MeasureType GetCurrentRegularizationEnergy();
	MeasureType GetCurrentEnergy();

	itkGetConstMacro( CurrentValue, MeasureType );

protected:
	/** Manual learning rate to apply. It is overridden by
	 * automatic learning rate estimation if enabled. See main documentation.
	 */
	InternalComputationValueType  m_LearningRate;

	/** Minimum convergence value for convergence checking.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the Functional profile. When the convergence value reaches
	 *  a small value, such as 1e-8, it would be treated as converged.
	 */
	InternalComputationValueType m_MinimumConvergenceValue;

	/** Window size for the convergence checker.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the Functional (metric value) profile.
	 */
	SizeValueType m_ConvergenceWindowSize;


	/** Current convergence value. */
	InternalComputationValueType m_ConvergenceValue;

	/** The convergence checker. */
	typename ConvergenceMonitoringType::Pointer m_ConvergenceMonitoring;

	SpectralOptimizer();
	~SpectralOptimizer() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	void SpectralUpdate( ParametersPointer parameters, const ParametersType* lambda, ParametersPointer nextParameters, bool changeDirection = false );

	virtual void InitializeAuxiliarParameters( void ) = 0;
	virtual void SetUpdate() = 0;
	virtual void Iterate(void) = 0;

//	/** Store the best value and related parameters */
//	MeasureType                  m_CurrentBestValue;
//	ParametersType               m_BestParameters;
//
//	/** Flag to control returning of best value and parameters. */
//	bool m_ReturnBestParametersAndValue;

	/** Particular parameter definitions from our method */
	InternalComputationValueType m_StepSize; // Step-size is tau in the formulations
	MatrixType m_A;
	MatrixType m_B;


	ParametersPointer m_Parameters;
	ParametersPointer m_NextParameters;
	TensorFieldPointer m_Denominator;
	MeasureType m_CurrentValue;

	FunctionalPointer m_Functional;
	MeasureType m_RegularizationEnergy;

	ParametersPointer m_LastField;
	//ParametersPointer m_CurrentField;

	/* Common variables for optimization control and reporting */
	bool                          m_Stop;
	StopConditionType             m_StopCondition;
	StopConditionDescriptionType  m_StopConditionDescription;
	SizeValueType                 m_NumberOfIterations;
	SizeValueType                 m_CurrentIteration;
	GridSizeType                  m_GridSize;

private:
	SpectralOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	void ApplyRegularizationTerm( ComplexFieldType* reference );
	void InitializeParameters( void );
	void InitializeDenominator( ComplexFieldType* reference );
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SpectralOptimizer.hxx"
#endif


#endif /* SPECTRALOPTIMIZER_H_ */
