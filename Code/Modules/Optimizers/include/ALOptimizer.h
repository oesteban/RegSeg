// --------------------------------------------------------------------------------------
// File:          ALOptimizer.h
// Date:          Jun 24, 2013
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

#ifndef ALOPTIMIZER_H_
#define ALOPTIMIZER_H_

#include <itkObject.h>
#include "LevelSetsOptimizer.h"
#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>
#include <itkImageIteratorWithIndex.h>
#include <itkImageAlgorithm.h>

namespace rstk
{
/**
 * \class ALOptimizer
 *  \brief Gradient descent optimizer.
 *
 * GradientDescentOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current deformation field is updated according:
 * \f[
 *        u^{t+1} = \mathcal{FT^{-1}}
 * \f]
 */

template< typename TLevelSetsFunction >
class ALOptimizer: public LevelSetsOptimizer<TLevelSetsFunction> {
public:
	/** Standard class typedefs and macros */
	typedef ALOptimizer          Self;
	typedef LevelSetsOptimizer<TLevelSetsFunction>     Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( ALOptimizer, LevelSetsOptimizer ); // Run-time type information (and related methods)
	itkNewMacro( Self );                                             // New macro for creation of through a Smart Pointer


	itkStaticConstMacro( Dimension, unsigned int, Superclass::Dimension );

	/** Metric type over which this class is templated */
	typedef typename Superclass::PointType                                           PointType;
	typedef typename Superclass::VectorType                                           VectorType;
	typedef typename Superclass::PointValueType                                      PointValueType;

	typedef typename Superclass::MeasureType                                         MeasureType;
	typedef typename Superclass::InternalComputationValueType                        InternalComputationValueType;
	typedef typename Superclass::StopConditionType                                   StopConditionType;
	typedef typename Superclass::LevelSetsFunctionType                               LevelSetsFunctionType;
	typedef typename Superclass::LevelSetsPointer                                    LevelSetsPointer;
	typedef typename Superclass::DeformationFieldType                                DeformationFieldType;
	typedef typename Superclass::DeformationFieldPointer                             DeformationFieldPointer;
	typedef typename DeformationFieldType::ConstPointer                              DeformationFieldConstPointer;
	typedef typename Superclass::DeformationFieldDirectionType                       DeformationFieldDirectionType;
	typedef typename Superclass::GridSizeType                                        GridSizeType;
	typedef typename Superclass::MatrixType                                          MatrixType;
	typedef typename Superclass::DeformationComponentType                            DeformationComponentType;
	typedef typename Superclass::DeformationComponentPointer                         DeformationComponentPointer;
	typedef typename Superclass::ContourType                                         ContourType;
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

	itkSetObjectMacro(uField, DeformationFieldType);
	itkGetConstObjectMacro(uField, DeformationFieldType);

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

	itkSetMacro( R, InternalComputationValueType );
	itkGetConstMacro( R, InternalComputationValueType );

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
	InternalComputationValueType m_R; // Step-size is tau in the formulations
	InternalComputationValueType m_Rho;
	MatrixType m_A;
	MatrixType m_B;


	DeformationFieldPointer m_uField;
	DeformationFieldPointer m_uFieldNext;
	DeformationFieldPointer m_vField;
	DeformationFieldPointer m_vFieldNext;
	DeformationFieldPointer m_lambdaField;
	DeformationFieldPointer m_lambdaFieldNext;

	DeformationFieldPointer m_NextDeformationField;
	TensorFieldPointer m_Denominator;
	MeasureType m_CurrentLevelSetsValue;


	ALOptimizer();
	~ALOptimizer() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	//void Iterate(void);
	void ApplyRegularizationTerm( ComplexFieldType* reference );
	void InitializeFields( const PointType orig, const PointType end, const DeformationFieldDirectionType dir);
	void InitializeDenominator( ComplexFieldType* reference );

	void UpdateU(void);
	void UpdateV(void);
	void UpdateLambda(void);
	void SetUpdate(void);
private:
	ALOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ALOptimizer.hxx"
#endif



#endif /* ALOPTIMIZER_H_ */
