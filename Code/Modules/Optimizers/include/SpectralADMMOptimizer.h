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

#ifndef SPECTRALADMMOPTIMIZER_H_
#define SPECTRALADMMOPTIMIZER_H_

#include "SpectralOptimizer.h"

using namespace itk;

namespace rstk
{
/**
 * \class SpectralADMMOptimizer
 * \brief Alternating Direction Method of Multipliers (ADMM) augmented
 * lagrangian descent optimizer.
 *
 * The alternating direction method of multipliers (ADMM) is a variant of the augmented
 * Lagrangian scheme that uses partial updates for the dual variables
 */

template< typename TFunctional >
class SpectralADMMOptimizer: public SpectralOptimizer<TFunctional> {
public:
	/** Standard class typedefs and macros */
	typedef SpectralADMMOptimizer                      Self;
	typedef SpectralOptimizer<TFunctional>             Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( SpectralADMMOptimizer, SpectralOptimizer ); // Run-time type information (and related methods)
	itkNewMacro( Self );                                             // New macro for creation of through a Smart Pointer

	/** Metric type over which this class is templated */
	typedef typename Superclass::FunctionalType                   FunctionalType;
	itkStaticConstMacro( Dimension, unsigned int, FunctionalType::Dimension );

	/** Codes of stopping conditions. */
	typedef typename Superclass::StopConditionType                StopConditionType;

	/** Stop condition return string type */
	typedef typename Superclass::StopConditionReturnStringType    StopConditionReturnStringType;

	/** Stop condition internal string type */
	typedef typename Superclass::StopConditionDescriptionType     StopConditionDescriptionType;

	/** Internal computation value type */
	typedef typename Superclass::InternalComputationValueType     InternalComputationValueType;

	/** Functional definitions */
	typedef typename Superclass::FunctionalPointer                FunctionalPointer;
	typedef typename Superclass::ParametersType                   ParametersType;
	typedef typename Superclass::DerivativeType                   DerivativeType;
	typedef typename Superclass::MeasureType                      MeasureType;
	typedef typename Superclass::PointType                        PointType;
	typedef typename Superclass::VectorType                       VectorType;
	typedef typename Superclass::PointValueType                   PointValueType;
	typedef typename Superclass::ParametersPointer                ParametersPointer;
	typedef typename Superclass::ParametersConstPointer           ParametersConstPointer;
	typedef typename Superclass::ParametersPointType              ParametersPointType;
	typedef typename Superclass::ParametersDirectionType          ParametersDirectionType;
	typedef typename Superclass::GridSizeType                     GridSizeType;
	typedef typename Superclass::ParametersComponentType          ParametersComponentType;
	typedef typename Superclass::ParametersComponentPointer       ParametersComponentPointer;
	typedef typename Superclass::FFTType                          FFTType;
	typedef typename Superclass::FFTPointer                       FFTPointer;
	typedef typename Superclass::FTDomainType                     FTDomainType;
	typedef typename Superclass::FTDomainPointer                  FTDomainPointer;
	typedef typename Superclass::ComplexType                      ComplexType;
	typedef typename Superclass::ComplexValueType                 ComplexValueType;
	typedef typename Superclass::IFFTType                         IFFTType;
	typedef typename Superclass::IFFTPointer                      IFFTPointer;
	typedef typename Superclass::RealPartType                     RealPartType;
	typedef typename Superclass::ComplexValuesVector              ComplexValuesVector;
	typedef typename Superclass::ComplexFieldValue                ComplexFieldValue;
	typedef typename Superclass::ComplexFieldType                 ComplexFieldType;
	typedef typename Superclass::ComplexFieldPointer              ComplexFieldPointer;
	typedef typename Superclass::MatrixType                       MatrixType;
	typedef typename Superclass::TensorFieldType                  TensorFieldType;
	typedef typename Superclass::TensorFieldPointer               TensorFieldPointer;
	typedef typename Superclass::SizeValueType                    SizeValueType;

	/** Type for the convergence checker */
	typedef typename Superclass::ConvergenceMonitoringType        ConvergenceMonitoringType;

	itkSetMacro( Rho, InternalComputationValueType);
	itkGetConstMacro( Rho, InternalComputationValueType);

protected:
	SpectralADMMOptimizer();
	~SpectralADMMOptimizer() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	void InitializeAuxiliarParameters( void );
	void Iterate(void);
	void SetUpdate();

private:
	SpectralADMMOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	void UpdateU(void);
	void UpdateV(void);
	void UpdateLambda(void);

	InternalComputationValueType m_Rho;
	ParametersPointer m_vField;
	ParametersPointer m_vFieldNext;
	ParametersPointer m_lambdaField;
	ParametersPointer m_lambdaFieldNext;

}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SpectralADMMOptimizer.hxx"
#endif


#endif /* SPECTRALADMMOPTIMIZER_H_ */
