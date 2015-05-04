// --------------------------------------------------------------------------------------
// File:          SegmentationOptimizer.h
// Date:          Feb 10, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of RegSeg
//
// RegSeg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RegSeg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RegSeg.  If not, see <http://www.gnu.org/licenses/>.
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
