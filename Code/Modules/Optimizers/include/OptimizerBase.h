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

#ifndef OPTIMIZERBASE_H_
#define OPTIMIZERBASE_H_


#include <boost/program_options.hpp>

#include <itkObjectToObjectOptimizerBase.h>

#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>

#include <itkImageIteratorWithIndex.h>
#include <itkImageAlgorithm.h>

#include "rstkMacro.h"
#include "ConfigurableObject.h"

using namespace itk;
namespace bpo = boost::program_options;


namespace rstk
{
/**
 * \class OptimizerBase
 *  \brief Gradient descent optimizer.
 *
 * GradientDescentOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current deformation field is updated according:
 * \f[
 *        u^{t+1} = \mathcal{FT^{-1}}
 * \f]
 */

template< typename TFunctional >
class OptimizerBase:
		public ObjectToObjectOptimizerBase,
		public ConfigurableObject
{
public:
	/** Standard class typedefs and macros */
	typedef OptimizerBase           Self;
	typedef ObjectToObjectOptimizerBase                Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;
	typedef ConfigurableObject                         SettingsClass;
	typedef typename SettingsClass::SettingsMap        SettingsMap;
	typedef typename SettingsClass::SettingsDesc       SettingsDesc;

	itkTypeMacro( OptimizerBase, ObjectToObjectOptimizerBase ); // Run-time type information (and related methods)

	/** Metric type over which this class is templated */
	typedef TFunctional                                             FunctionalType;
	typedef typename FunctionalType::ScalesType                     GradientScales;
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

	typedef size_t                                                  SizeValueType;
	typedef typename Superclass::ScalesType                         ScalesType;

	/** Functional definitions */
	typedef typename FunctionalType::Pointer                        FunctionalPointer;
	typedef typename FunctionalType::MeasureType                    MeasureType;
	typedef itk::Array< MeasureType >                               MeasureArray;
	typedef typename FunctionalType::PointType                      PointType;
	typedef typename FunctionalType::VectorType                     VectorType;
	typedef typename FunctionalType::PointValueType                 PointValueType;
	typedef typename FunctionalType::PointValueType                 InternalComputationValueType; // This comes from MetricType in Superclass


	typedef SparseMatrixTransform<PointValueType, Dimension >       TransformType;
	typedef typename TransformType::Pointer                         TransformPointer;
	typedef typename TransformType::CoefficientsImageType           CoefficientsImageType;
	typedef typename TransformType::CoeffImagePointer               CoefficientsImagePointer;
	typedef typename TransformType::CoefficientsImageArray          CoefficientsImageArray;
	typedef typename TransformType::ParametersType                  ParametersType;
	typedef typename TransformType::WeightsMatrix                   WeightsMatrix;
	typedef typename TransformType::DimensionVector					ParametersVector;
	typedef typename TransformType::DimensionParameters             ParametersContainer;
	typedef itk::FixedArray< ParametersVector*, Dimension >         ParametersPointerContainer;
	typedef typename TransformType::FieldType                       FieldType;
	typedef typename TransformType::FieldPointer                    FieldPointer;
	typedef typename TransformType::FieldConstPointer               FieldConstPointer;
	typedef typename TransformType::SizeType                        ControlPointsGridSizeType;
	typedef typename TransformType::SpacingType                     ControlPointsGridSpacingType;

	/** Type for the convergence checker */
	typedef itk::Function::
			WindowConvergenceMonitoringFunction<double>	            ConvergenceMonitoringType;

	/** Accessors for Functional */
	itkGetObjectMacro( Functional, FunctionalType );
	itkSetObjectMacro( Functional, FunctionalType );

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

	/** Get stop condition enum */
	itkGetConstReferenceMacro(StopCondition, StopConditionType);

	itkGetConstMacro(CurrentValue, MeasureType);
	itkGetConstMacro(CurrentNorm, MeasureType);

	itkSetMacro( DescriptorRecompPeriod, SizeValueType );
	itkGetConstMacro( DescriptorRecompPeriod, SizeValueType );

	itkSetMacro( UseDescriptorRecomputation, bool );
	itkGetConstMacro( UseDescriptorRecomputation, bool );

	itkGetConstMacro( GridSize, ControlPointsGridSizeType );
	itkGetConstMacro( GridSpacing, ControlPointsGridSpacingType );

	itkGetConstMacro( MaxDisplacement, VectorType );

	itkGetConstMacro( IsDiffeomorphic, bool );
	itkGetConstMacro( DiffeomorphismForced, bool );

	itkGetConstMacro( ForceDiffeomorphic, bool );
	itkSetMacro( ForceDiffeomorphic, bool );

	itkSetMacro(LearningRate, InternalComputationValueType);               // Set the learning rate
	itkGetConstReferenceMacro(LearningRate, InternalComputationValueType); // Get the learning rate

	itkSetMacro( AutoStepSize, bool );
	itkGetConstMacro( AutoStepSize, bool );

	itkSetMacro( UseLightWeightConvergenceChecking, bool );
	itkGetConstMacro( UseLightWeightConvergenceChecking, bool );

	virtual void SetStepSize (const InternalComputationValueType _arg);
	itkGetConstMacro( StepSize, InternalComputationValueType );
	itkGetConstMacro( Momentum, InternalComputationValueType );
	itkGetConstMacro( MaximumGradient, InternalComputationValueType );

	itkGetConstMacro( MaxSpeed, InternalComputationValueType );
	itkGetConstMacro( AvgSpeed, InternalComputationValueType );
	itkGetConstMacro( MeanSpeed, InternalComputationValueType );

	itkSetMacro(Coefficients, CoefficientsImageArray);
	itkGetConstMacro(Coefficients, CoefficientsImageArray);

	itkSetMacro(DerivativeCoefficients, CoefficientsImageArray);
	itkGetConstMacro(DerivativeCoefficients, CoefficientsImageArray);

	itkGetObjectMacro( Transform, TransformType );

	virtual const FieldType * GetCurrentCoefficients() const = 0;
	virtual const FieldType * GetCurrentCoefficientsField() const = 0;

	/** Start and run the optimization */
	void Start();

	void Stop(void);

	/** Get the reason for termination */
	const StopConditionReturnStringType GetStopConditionDescription() const;

	void Resume();

	virtual MeasureType GetCurrentRegularizationEnergy() = 0;
	virtual MeasureType GetCurrentEnergy() = 0;

	static void AddOptions( SettingsDesc& opts );
protected:
	OptimizerBase();
	~OptimizerBase() {}
	void PrintSelf( std::ostream &os, itk::Indent indent ) const;
	virtual void ParseSettings();

	virtual void InitializeParameters() = 0;
	virtual void InitializeAuxiliarParameters() = 0;
	virtual void ComputeDerivative() = 0;
	virtual void Iterate() = 0;
	virtual void PostIteration() = 0;

	virtual bool DoDescriptorsUpdate();

	/** Manual learning rate to apply. It is overridden by
	 * automatic learning rate estimation if enabled. See main documentation.
	 */
	InternalComputationValueType  m_LearningRate;
	InternalComputationValueType  m_Momentum;
	InternalComputationValueType  m_LastMaximumGradient;
	InternalComputationValueType  m_MaximumGradient;

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


	/* Common variables for optimization control and reporting */
	bool                          m_Stop;
	StopConditionType             m_StopCondition;
	StopConditionDescriptionType  m_StopConditionDescription;
	SizeValueType                 m_DescriptorRecompPeriod;
	SizeValueType                 m_NextRecompIteration;
	SizeValueType                 m_ValueOscillations;
	SizeValueType                 m_ValueOscillationsMax;
	SizeValueType                 m_ValueOscillationsLast;
	bool                          m_UseDescriptorRecomputation;

	InternalComputationValueType  m_StepSize;
	InternalComputationValueType  m_MaxSpeed;
	InternalComputationValueType  m_MeanSpeed;
	InternalComputationValueType  m_AvgSpeed;
	bool                          m_AutoStepSize;
	bool                          m_IsDiffeomorphic;
	bool                          m_ForceDiffeomorphic;
	bool                          m_DiffeomorphismForced;
	bool                          m_UseLightWeightConvergenceChecking;
	bool                          m_UseAdaptativeDescriptors;

	/* Energy tracking */
	MeasureType                  m_CurrentValue;
	MeasureType                  m_CurrentEnergy;
	MeasureType                  m_CurrentNorm;
	MeasureType                  m_LastEnergy;
	MeasureArray                 m_ValueWindow;

	ControlPointsGridSizeType    m_GridSize;
	ControlPointsGridSpacingType m_GridSpacing;
	VectorType                   m_MaxDisplacement;

	TransformPointer             m_Transform;
	FunctionalPointer            m_Functional;

	CoefficientsImageArray       m_Coefficients;
	CoefficientsImageArray       m_DerivativeCoefficients;
private:
	OptimizerBase( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class


itkEventMacro(FunctionalModifiedEvent, itk::AnyEvent);

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "OptimizerBase.hxx"
#endif




#endif /* OPTIMIZERBASE_H_ */
