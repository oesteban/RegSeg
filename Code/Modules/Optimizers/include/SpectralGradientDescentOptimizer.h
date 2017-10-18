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

#ifndef SPECTRALGRADIENTDESCENTOPTIMIZER_H_
#define SPECTRALGRADIENTDESCENTOPTIMIZER_H_

#include "SpectralOptimizer.h"

using namespace itk;

namespace rstk
{
/**
 * \class SpectralGradientDescentOptimizer
 *  \brief Gradient descent optimizer.
 *
 * GradientDescentOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current deformation field is updated according:
 * \f[
 *        u^{t+1} = \mathcal{FT^{-1}}
 * \f]
 */

template< typename TFunctional >
class SpectralGradientDescentOptimizer: public SpectralOptimizer<TFunctional> {
public:
	/** Standard class typedefs and macros */
	typedef SpectralGradientDescentOptimizer           Self;
	typedef SpectralOptimizer<TFunctional>             Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	itkTypeMacro( SpectralGradientDescentOptimizer, SpectralOptimizer ); // Run-time type information (and related methods)
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

	/** Functional definitions */
	typedef typename Superclass::FunctionalPointer                FunctionalPointer;
	typedef typename Superclass::ParametersType                   ParametersType;
	typedef typename Superclass::MeasureType                      MeasureType;
	typedef typename Superclass::PointType                        PointType;
	typedef typename Superclass::VectorType                       VectorType;
	typedef typename Superclass::PointValueType                   PointValueType;
	typedef typename Superclass::FFTType                          FFTType;
	typedef typename Superclass::FFTPointer                       FFTPointer;
	typedef typename Superclass::FTDomainType                     FTDomainType;
	typedef typename Superclass::FTDomainPointer                  FTDomainPointer;
	typedef typename Superclass::ComplexType                      ComplexType;
	typedef typename Superclass::InternalComputationValueType     InternalComputationValueType;
	typedef typename Superclass::InternalVectorType               InternalVectorType;
	typedef typename Superclass::InternalVectorFieldType          InternalVectorFieldType;
	typedef typename Superclass::InternalVectorFieldPointer       InternalVectorFieldPointer;

	typedef typename Superclass::IFFTType                         IFFTType;
	typedef typename Superclass::IFFTPointer                      IFFTPointer;
	typedef typename Superclass::RealPartType                     RealPartType;
	typedef typename Superclass::ComplexFieldValue                ComplexFieldValue;
	typedef typename Superclass::ComplexFieldType                 ComplexFieldType;
	typedef typename Superclass::ComplexFieldPointer              ComplexFieldPointer;
//	typedef typename Superclass::MatrixType                       MatrixType;
//	typedef typename Superclass::TensorFieldType                  TensorFieldType;
//	typedef typename Superclass::TensorFieldPointer               TensorFieldPointer;
	typedef typename Superclass::SizeValueType                    SizeValueType;

	/** Type for the convergence checker */
	typedef typename Superclass::ConvergenceMonitoringType        ConvergenceMonitoringType;

	typedef typename Superclass::CoefficientsImageType            CoefficientsImageType;
protected:
	SpectralGradientDescentOptimizer();
	~SpectralGradientDescentOptimizer() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	void InitializeAuxiliarParameters( void ) {}
	void Iterate(void);
	void SetUpdate();

private:
	SpectralGradientDescentOptimizer( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SpectralGradientDescentOptimizer.hxx"
#endif

#endif /* SPECTRALGRADIENTDESCENTOPTIMIZER_H_ */
