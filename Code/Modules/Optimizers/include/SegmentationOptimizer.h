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

#ifndef SEGMENTATIONOPTIMIZER_H_
#define SEGMENTATIONOPTIMIZER_H_


#include <boost/program_options.hpp>

#include <itkObjectToObjectOptimizerBase.h>

#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkRealToHalfHermitianForwardFFTImageFilter.h>
#include <itkHalfHermitianToRealInverseFFTImageFilter.h>

#include <itkImageIteratorWithIndex.h>
#include <itkImageAlgorithm.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>

#include "rstkMacro.h"
#include "ConfigurableObject.h"
#include "BSplineSparseMatrixTransform.h"

using namespace itk;
namespace bpo = boost::program_options;


namespace rstk
{
/**
 * \class SegmentationOptimizer
 *  \brief Gradient descent optimizer.
 *
 * GradientDescentOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current deformation field is updated according:
 * \f[
 *        u^{t+1} = \mathcal{FT^{-1}}
 * \f]
 */

template< typename TFunctional >
class SegmentationOptimizer: public OptimizerBase< TFunctional >
{
public:
	/** Standard class typedefs and macros */
	typedef SegmentationOptimizer           Self;
	typedef OptimizerBase< TFunctional >               Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	typedef typename Superclass::SettingsClass         SettingsClass;
	typedef typename Superclass::SettingsMap           SettingsMap;
	typedef typename Superclass::SettingsDesc          SettingsDesc;

	itkTypeMacro( SegmentationOptimizer, OptimizerBase ); // Run-time type information (and related methods)
	itkNewMacro( Self );                                  // New macro for creation of through a Smart Pointer

	/** Metric type over which this class is templated */
	typedef typename Superclass::FunctionalType        FunctionalType;
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
	typedef typename Superclass::StopConditionReturnStringType      StopConditionReturnStringType;

	/** Stop condition internal string type */
	typedef typename Superclass::StopConditionDescriptionType       StopConditionDescriptionType;
	typedef typename Superclass::SizeValueType                      SizeValueType;

	/** Functional definitions */
	typedef typename FunctionalType::SampleType                     SampleType;
	typedef typename FunctionalType::GradientSample                 GradientSample;
	typedef typename Superclass::FunctionalPointer                  FunctionalPointer;
	typedef typename Superclass::MeasureType                        MeasureType;
	typedef typename Superclass::PointType                          PointType;
	typedef typename Superclass::VectorType                         VectorType;
	typedef typename Superclass::PointValueType                     PointValueType;
	typedef typename Superclass::TransformType                      TransformType;
	typedef typename Superclass::TransformPointer                   TransformPointer;
	typedef typename Superclass::CoefficientsImageType              CoefficientsImageType;
	typedef typename Superclass::CoefficientsImageArray             CoefficientsImageArray;
	typedef typename Superclass::ParametersType                     ParametersType;
	typedef typename Superclass::WeightsMatrix                      WeightsMatrix;
	typedef typename Superclass::FieldType                          FieldType;
	typedef typename Superclass::FieldPointer                       FieldPointer;
	typedef typename Superclass::FieldConstPointer                  FieldConstPointer;

	/** Type for the convergence checker */
	typedef typename Superclass::ConvergenceMonitoringType          ConvergenceMonitoringType;

	MeasureType GetCurrentRegularizationEnergy() { return 0.0; };
	MeasureType GetCurrentEnergy() { return 0.0; };

	const FieldType * GetCurrentDisplacementField () const {
		return FieldType::New();
	}

protected:
	SegmentationOptimizer();
	~SegmentationOptimizer(){}
	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	void ParseSettings() { Superclass::ParseSettings(); }

	void InitializeParameters();
	void InitializeAuxiliarParameters() {};
	void ComputeDerivative();
	void Iterate();
	void PostIteration();
private:
	SegmentationOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	WeightsMatrix m_Gradient;
	WeightsMatrix m_Displacement;
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SegmentationOptimizer.hxx"
#endif




#endif /* SEGMENTATIONOPTIMIZER_H_ */
