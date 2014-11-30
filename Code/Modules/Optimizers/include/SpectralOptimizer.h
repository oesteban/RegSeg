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

#include <boost/program_options.hpp>

#include <itkWindowConvergenceMonitoringFunction.h>
#include <vector>
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkRealToHalfHermitianForwardFFTImageFilter.h>
#include <itkHalfHermitianToRealInverseFFTImageFilter.h>

#include <itkImageIteratorWithIndex.h>
#include <itkImageAlgorithm.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkAddImageFilter.h>


#include "rstkMacro.h"
#include "OptimizerBase.h"
#include "BSplineSparseMatrixTransform.h"

using namespace itk;
namespace bpo = boost::program_options;


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
class SpectralOptimizer: public OptimizerBase< TFunctional >
{
public:
	/** Standard class typedefs and macros */
	typedef SpectralOptimizer                          Self;
	typedef OptimizerBase< TFunctional >               Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( SpectralOptimizer, OptimizerBase ); // Run-time type information (and related methods)

	/* Configurable object typedefs */
	typedef typename Superclass::SettingsClass         SettingsClass;
	typedef typename Superclass::SettingsMap           SettingsMap;
	typedef typename Superclass::SettingsDesc          SettingsDesc;

	/** Metric type over which this class is templated */
	typedef TFunctional                                FunctionalType;
	itkStaticConstMacro( Dimension, unsigned int, FunctionalType::Dimension );

	/** Codes of stopping conditions. */
	using Superclass::StopConditionType;

	/** Inherited definitions */
	typedef typename Superclass::StopConditionReturnStringType  StopConditionReturnStringType;
	typedef typename Superclass::StopConditionDescriptionType   StopConditionDescriptionType;
    typedef typename Superclass::SizeValueType                  SizeValueType;
	typedef typename Superclass::ConvergenceMonitoringType	    ConvergenceMonitoringType;

	typedef typename Superclass::FunctionalPointer              FunctionalPointer;
	typedef typename Superclass::MeasureType                    MeasureType;
	typedef typename Superclass::PointType                      PointType;
	typedef typename Superclass::VectorType                     VectorType;
	typedef typename Superclass::PointValueType                 PointValueType;

	typedef typename Superclass::TransformType                  TransformType;
	typedef typename Superclass::TransformPointer               TransformPointer;
	typedef typename Superclass::CoefficientsImageType          CoefficientsImageType;
	typedef typename CoefficientsImageType::PixelType           CoefficientsValueType;
	typedef typename Superclass::CoefficientsImagePointer       CoefficientsImagePointer;
	typedef typename Superclass::CoefficientsImageArray         CoefficientsImageArray;
	typedef typename Superclass::ParametersType                 ParametersType;
	typedef typename Superclass::WeightsMatrix                  WeightsMatrix;
	typedef typename Superclass::ParametersVector				ParametersVector;
	typedef typename Superclass::ParametersContainer            ParametersContainer;
	typedef typename Superclass::FieldType                      FieldType;
	typedef typename Superclass::FieldPointer                   FieldPointer;
	typedef typename Superclass::FieldConstPointer              FieldConstPointer;
	typedef typename Superclass::ControlPointsGridSizeType      ControlPointsGridSizeType;
	typedef typename Superclass::ControlPointsGridSpacingType   ControlPointsGridSpacingType;

	typedef itk::MultiplyImageFilter<CoefficientsImageType, CoefficientsImageType, CoefficientsImageType> MultiplyFilterType;
	typedef itk::AddImageFilter<CoefficientsImageType, CoefficientsImageType, CoefficientsImageType>      AddFilterType;
	typedef itk::AddImageFilter<FieldType, FieldType, FieldType> AddFieldFilterType;
	typedef typename AddFieldFilterType::Pointer                 AddFieldFilterPointer;

	typedef BSplineSparseMatrixTransform
			                      < PointValueType, Dimension, 3u > SplineTransformType;
	typedef typename SplineTransformType::Pointer                   SplineTransformPointer;

	typedef itk::RealToHalfHermitianForwardFFTImageFilter
			                          <CoefficientsImageType>       FFTType;
	typedef typename FFTType::Pointer                               FFTPointer;
	typedef typename FFTType::OutputImageType                       FTDomainType;
	typedef typename FTDomainType::Pointer                          FTDomainPointer;
	typedef itk::FixedArray< FTDomainPointer, Dimension >           FTDomainArray;
	typedef typename FTDomainType::PixelType                        ComplexType;
	typedef itk::AddImageFilter<FTDomainType, FTDomainType, FTDomainType>
																	FTAddFilterType;
	typedef itk::MultiplyImageFilter<FTDomainType, FTDomainType, FTDomainType>
																	FTMultiplyFilterType;
	typedef itk::DivideImageFilter<FTDomainType, FTDomainType, FTDomainType>
																	FTDivideFilterType;

	/** Internal computation value type */
	typedef typename ComplexType::value_type                        InternalComputationValueType;
	typedef itk::Vector< InternalComputationValueType, Dimension >  InternalVectorType;
	typedef itk::Image< InternalVectorType, Dimension >             InternalVectorFieldType;
	typedef typename InternalVectorFieldType::Pointer               InternalVectorFieldPointer;
	typedef itk::ContinuousIndex< InternalComputationValueType, Dimension>
																	ContinuousIndexType;

	typedef itk::HalfHermitianToRealInverseFFTImageFilter
			        <FTDomainType, CoefficientsImageType>           IFFTType;
	typedef typename IFFTType::Pointer                              IFFTPointer;


	typedef itk::Image< InternalComputationValueType, Dimension >   RealPartType;
	typedef itk::Vector< ComplexType, Dimension >                   ComplexFieldValue;
	typedef itk::Image< ComplexFieldValue, Dimension >              ComplexFieldType;
	typedef typename ComplexFieldType::Pointer                      ComplexFieldPointer;

	itkSetMacro( Alpha, InternalVectorType );
	itkGetConstMacro( Alpha, InternalVectorType );

	itkSetMacro( Beta, InternalVectorType );
	itkGetConstMacro( Beta, InternalVectorType );

	void SetAlpha(const InternalComputationValueType v ) {
		this->m_Alpha.Fill(v);
		this->Modified();
	}

	void SetBeta(const InternalComputationValueType v ) {
		this->m_Beta.Fill(v);
		this->Modified();
	}

	itkSetMacro( GridSize, ControlPointsGridSizeType );
	// void SetGridSize( double val ) { this->m_GridSize.Fill(val); }

	itkSetMacro( GridSpacing, ControlPointsGridSpacingType );
	// void SetGridSpacing( double val ) { this->m_GridSize.Fill(val); }

	void ComputeIterationSpeed();
	MeasureType GetCurrentRegularizationEnergy();
	MeasureType GetCurrentEnergy();

	itkSetMacro( DescriptorRecomputationFreq, SizeValueType );
	itkGetConstMacro( DescriptorRecomputationFreq, SizeValueType );

	itkSetMacro( UseDescriptorRecomputation, bool );
	itkGetConstMacro( UseDescriptorRecomputation, bool );

	itkGetConstObjectMacro(CurrentCoefficients, FieldType);
	itkSetObjectMacro(InitialDisplacementField, FieldType);

	const FieldType * GetCurrentCoefficientsField () const {
		return this->m_Transform->GetCoefficientsVectorImage();
	}

	static void AddOptions( SettingsDesc& opts );
protected:
	SpectralOptimizer();
	~SpectralOptimizer() {}
	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	/* Inherited from OptimizerBase */
	virtual void ComputeDerivative();
	virtual void Iterate() = 0;
	virtual void PostIteration();
	void InitializeParameters();
	virtual void InitializeAuxiliarParameters() = 0;

	/* SpectralOptimizer specific members */
	void ComputeUpdate(CoefficientsImageArray uk,
			           const CoefficientsImageArray gk,
			           CoefficientsImageArray next_uk,
			           bool changeDirection = false );
	void BetaRegularization(CoefficientsImagePointer numerator,
			                CoefficientsImageArray next_uk,
			                InternalComputationValueType s,
			                size_t d);

	virtual void SetUpdate() = 0;

	virtual void ParseSettings();

	/* Common variables for optimization control and reporting */
	bool                          m_DenominatorCached;

	/** Particular parameter definitions from our method */
	InternalVectorType m_Alpha;
	InternalVectorType m_Beta;

	/* Energy tracking */
	MeasureType m_RegularizationEnergy;
	MeasureType m_CurrentTotalEnergy;
	bool m_RegularizationEnergyUpdated;

	CoefficientsImageArray       m_NextCoefficients;
	CoefficientsImageArray       m_Denominator;
	FTDomainPointer              m_FTLaplacian;
	FTDomainPointer              m_FTOnes;
	FieldPointer                 m_LastCoeff;
	FieldPointer                 m_InitialCoeff;
	FieldPointer                 m_CurrentCoefficients;
	FieldPointer                 m_InitialDisplacementField;
	AddFieldFilterPointer        m_FieldCoeffAdder;
private:
	SpectralOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	void ApplyRegularizationTerm( ComplexFieldType* reference );
	void ApplyRegularizationComponent( size_t d, FTDomainType *reference );
	void InitializeDenominator( itk::ImageBase<Dimension> *reference );
	void UpdateField();
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SpectralOptimizer.hxx"
#endif


#endif /* SPECTRALOPTIMIZER_H_ */
