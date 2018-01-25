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

	using typename SparseMatrixTransform< TScalar, NDimensions >::InterpolateModeType;
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
