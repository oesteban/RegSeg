// --------------------------------------------------------------------------------------
// File:          SpectralADMMOptimizer.h
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
