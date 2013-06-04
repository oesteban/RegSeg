/* --------------------------------------------------------------------------------------
 * File:    GradientDescentLevelSetsOptimizer.h
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

#ifndef GradientDescentLevelSetsOptimizer_h
#define GradientDescentLevelSetsOptimizer_h

#include <itkObject.h>
#include "LevelSetsOptimizer.h"
#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>
#include <itkImageIteratorWithIndex.h>
#include <itkImageAlgorithm.h>

namespace rstk
{
/**
 * \class GradientDescentLevelSetsOptimizer
 *  \brief Gradient descent optimizer.
 *
 * GradientDescentOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current deformation field is updated according:
 * \f[
 *        u^{t+1} = \mathcal{FT^{-1}}
 * \f]
 */

template< typename TLevelSetsFunction >
class GradientDescentLevelSetsOptimizer: public LevelSetsOptimizer<TLevelSetsFunction> {
public:
	/** Standard class typedefs and macros */
	typedef GradientDescentLevelSetsOptimizer          Self;
	typedef LevelSetsOptimizer<TLevelSetsFunction>     Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( GradientDescentLevelSetsOptimizer, LevelSetsOptimizer ); // Run-time type information (and related methods)
	itkNewMacro( Self );                                             // New macro for creation of through a Smart Pointer


	itkStaticConstMacro( Dimension, unsigned int, Superclass::Dimension );

	/** Metric type over which this class is templated */
	typedef typename Superclass::MeasureType                                         MeasureType;
	typedef typename Superclass::InternalComputationValueType                        InternalComputationValueType;
	typedef typename Superclass::StopConditionType                                   StopConditionType;
	typedef typename Superclass::LevelSetsFunctionType                               LevelSetsFunctionType;
	typedef typename Superclass::LevelSetsPointer                                    LevelSetsPointer;
	typedef typename Superclass::DeformationFieldType                                DeformationFieldType;
	typedef typename Superclass::DeformationFieldPointer                             DeformationFieldPointer;
	typedef typename Superclass::DeformationFieldPointType                           DeformationFieldPointType;
	typedef typename Superclass::DeformationFieldDirectionType                       DeformationFieldDirectionType;
	typedef typename Superclass::PointValueType                                      PointValueType;
	typedef typename Superclass::VectorType                                          VectorType;
	typedef typename Superclass::MatrixType                                          MatrixType;
	typedef typename Superclass::DeformationComponentType                            DeformationComponentType;
	typedef typename Superclass::DeformationComponentPointer                         DeformationComponentPointer;
	typedef typename Superclass::ContourDeformationType                              ContourDeformationType;
	typedef typename Superclass::ContourDeformationPointer                           ContourDeformationPointer;
	typedef typename Superclass::ContourPointType                                    ContourPointType;
	typedef typename Superclass::TensorFieldType                                     TensorFieldType;
	typedef typename Superclass::TensorFieldPointer                                  TensorFieldPointer;

	typedef typename Superclass::FFTType                                             FFTType;
	typedef typename Superclass::FFTPointer                                          FFTPointer;
	typedef typename Superclass::FTDomainType                                        FTDomainType;
	typedef typename Superclass::FTDomainPointer                                     FTDomainPointer;
	typedef typename Superclass::ComplexType                                         ComplexType;
	typedef typename Superclass::ComplexValueType                                    ComplexValueType;
	typedef typename Superclass::RealPartType                                        RealPartType;
	typedef typename RealPartType::Pointer                                           RealPartPointer;

	typedef typename Superclass::ComplexValuesVector                                 ComplexValuesVector;
	typedef typename Superclass::ComplexFieldValue                                   ComplexFieldValue;
	typedef typename Superclass::ComplexFieldType                                    ComplexFieldType;
	typedef typename Superclass::ComplexFieldPointer                                 ComplexFieldPointer;

	typedef typename Superclass::IFFTType                                            IFFTType;
	typedef typename Superclass::IFFTPointer                                         IFFTPointer;

	/** Type for the convergence checker */
	typedef itk::Function::WindowConvergenceMonitoringFunction<MeasureType>	         ConvergenceMonitoringType;
	typedef typename ConvergenceMonitoringType::EnergyValueContainerSizeType         SizeValueType;

	itkSetMacro(LearningRate, InternalComputationValueType);               // Set the learning rate
	itkGetConstReferenceMacro(LearningRate, InternalComputationValueType); // Get the learning rate

	itkSetObjectMacro(DeformationField, DeformationFieldType);
	itkGetConstObjectMacro(DeformationField, DeformationFieldType);

	/** Minimum convergence value for convergence checking.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the LevelSets profile. When the convergence value reaches
	 *  a small value, it would be treated as converged.
	 *
	 *  The default m_MinimumConvergenceValue is set to 1e-8 to pass all
	 *  tests. It is suggested to use 1e-6 for less stringent convergence
	 *  checking.
	 */
	itkSetMacro(MinimumConvergenceValue, InternalComputationValueType);

	/** Window size for the convergence checker.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the LevelSets (metric value) profile.
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

	/** Start and run the optimization */
	virtual void Start();

	virtual void Stop(void);

	virtual void Resume();

	virtual void ComputeIterationChange();

protected:
	/** Manual learning rate to apply. It is overridden by
	 * automatic learning rate estimation if enabled. See main documentation.
	 */
	InternalComputationValueType  m_LearningRate;

	/** The maximum step size in physical units, to restrict learning rates.
	 * Only used with automatic learning rate estimation. See main documentation.
	 */
	InternalComputationValueType  m_MaximumStepSizeInPhysicalUnits;
	/** Minimum convergence value for convergence checking.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the LevelSets profile. When the convergence value reaches
	 *  a small value, such as 1e-8, it would be treated as converged.
	 */
	InternalComputationValueType m_MinimumConvergenceValue;

	/** Window size for the convergence checker.
	 *  The convergence checker calculates convergence value by fitting to
	 *  a window of the LevelSets (metric value) profile.
	 */
	SizeValueType m_ConvergenceWindowSize;


	/** Current convergence value. */
	InternalComputationValueType m_ConvergenceValue;

	/** The convergence checker. */
	typename ConvergenceMonitoringType::Pointer m_ConvergenceMonitoring;

//	/** Store the best value and related paramters */
//	MeasureType                  m_CurrentBestValue;
//	ParametersType               m_BestParameters;
//
//	/** Flag to control returning of best value and parameters. */
//	bool m_ReturnBestParametersAndValue;

	/** Particular parameter definitions from our method */
	InternalComputationValueType m_StepSize; // Step-size is tau in the formulations
	MatrixType m_A;
	MatrixType m_B;


	DeformationFieldPointer m_DeformationField;
	DeformationFieldPointer m_NextDeformationField;
	DeformationFieldPointer m_ShapeGradients;
	TensorFieldPointer m_Denominator;
	MeasureType m_CurrentLevelSetsValue;


	GradientDescentLevelSetsOptimizer();
	~GradientDescentLevelSetsOptimizer() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	void Iterate(void);
	void ApplyRegularizationTerm( ComplexFieldType* reference );
	void InitializeDeformationField( const DeformationFieldPointType orig, const DeformationFieldPointType end, const DeformationFieldDirectionType dir);
	void InitializeDenominator( ComplexFieldType* reference );

private:
	GradientDescentLevelSetsOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "GradientDescentLevelSetsOptimizer.hxx"
#endif

#endif /* GradientDescentLevelSetsOptimizer_h */
