// --------------------------------------------------------------------------------------
// File:          SpectralGradientDescentOptimizer.h
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
