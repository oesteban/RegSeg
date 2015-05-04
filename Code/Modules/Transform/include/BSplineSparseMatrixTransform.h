// --------------------------------------------------------------------------------------
// File:          BSplineSparseMatrixTransform.h
// Date:          Jan 15, 2014
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

#ifndef BSPLINESPARSEMATRIXTRANSFORM_H_
#define BSPLINESPARSEMATRIXTRANSFORM_H_

#include "SparseMatrixTransform.h"
#include <itkBSplineKernelFunction.h>
#include <itkBSplineDerivativeKernelFunction.h>
#include "BSplineSecondDerivativeKernelFunction.h"


namespace rstk {

template< class TScalar, unsigned int NDimensions = 3u, unsigned int VSplineOrder = 3u >
class BSplineSparseMatrixTransform: public SparseMatrixTransform< TScalar, NDimensions > {
public:
	typedef BSplineSparseMatrixTransform< TScalar, NDimensions, VSplineOrder >          Self;
	typedef SparseMatrixTransform< TScalar, NDimensions >                               Superclass;
	typedef itk::SmartPointer< Self >                                                   Pointer;
	typedef itk::SmartPointer< const Self >                                             ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( BSplineSparseMatrixTransform, SparseMatrixTransform );
	itkStaticConstMacro( Dimension, unsigned int, NDimensions );
	itkStaticConstMacro( SplineOrder, unsigned int, VSplineOrder );

	typedef typename Superclass::ScalarType                                             ScalarType;
	typedef typename Superclass::PointType                                              PointType;
	typedef typename Superclass::VectorType                                             VectorType;
	typedef typename Superclass::KernelFunctionType                                     KernelFunctionType;
	typedef typename Superclass::KernelFunctionPointer                                  KernelFunctionPointer;
	typedef typename Superclass::WeightsMatrix                                          WeightsMatrix;
	typedef typename Superclass::DimensionVector                                        DimensionVector;
	typedef typename Superclass::ParametersType                                         ParametersType;
	typedef typename Superclass::PointsList                                             PointsList;
	typedef typename Superclass::JacobianType                                           JacobianType;
	typedef typename Superclass::CoefficientsImageArray                                 CoefficientsImageArray;
	typedef typename Superclass::CoefficientsImageType                                  CoefficientsImageType;

	/** Standard coordinate point type for this class. */
	typedef typename Superclass::InputPointType                                         InputPointType;
	typedef typename Superclass::OutputPointType                                        OutputPointType;

	/** Standard vector type for this class. */
	typedef typename Superclass::InputVectorType                                        InputVectorType;
	typedef typename Superclass::OutputVectorType                                       OutputVectorType;

	/** Standard covariant vector type for this class */
	typedef typename Superclass::InputCovariantVectorType                               InputCovariantVectorType;
	typedef typename Superclass::OutputCovariantVectorType                              OutputCovariantVectorType;

	/** Standard vnl_vector type for this class. */
	typedef typename Superclass::InputVnlVectorType                                     InputVnlVectorType;
	typedef typename Superclass::OutputVnlVectorType                                    OutputVnlVectorType;

	typedef typename Superclass::ArrayType                                              ArrayType;

	typedef typename Superclass::AltCoeffType                        AltCoeffType;
	typedef typename Superclass::AltCoeffPointer                     AltCoeffPointer;

	using Superclass::InterpolateModeType;
protected:
	BSplineSparseMatrixTransform(): Superclass() {
		this->m_KernelFunction = dynamic_cast< KernelFunctionType * >(
				itk::BSplineKernelFunction<SplineOrder, ScalarType>::New().GetPointer() );
		this->m_DerivativeKernel = dynamic_cast< KernelFunctionType * >(
				itk::BSplineDerivativeKernelFunction<SplineOrder, ScalarType>::New().GetPointer() );
		//this->m_SecondDerivativeKernel = dynamic_cast< KernelFunctionType * >(
		//        itk::BSplineSecondDerivativeKernelFunction<SplineOrder, ScalarType>::New().GetPointer() );
	}

	~BSplineSparseMatrixTransform() {}

	inline size_t GetSupport() const {
		return SplineOrder;
	}

private:
	BSplineSparseMatrixTransform( const Self & );
	void operator=( const Self & );
};

} // namespace rstk

#endif /* BSPLINESPARSEMATRIXTRANSFORM_H_ */
