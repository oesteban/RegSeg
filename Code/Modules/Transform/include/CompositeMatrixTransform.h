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

#ifndef COMPOSITEMATRIXTRANSFORM_H_
#define COMPOSITEMATRIXTRANSFORM_H_

#include <vector>
#include <iostream>
#include <itkDisplacementFieldTransform.h>
#include <itkMatrix.h>

#include "CachedMatrixTransform.h"
#include "SparseMatrixTransform.h"
#include "BSplineSparseMatrixTransform.h"
#include "rstkMacro.h"

namespace rstk {

template< class TScalar, unsigned int NDimensions = 3u >
class CompositeMatrixTransform : public rstk::CachedMatrixTransform< TScalar, NDimensions >
{
public:
    /* Standard class typedefs. */
    typedef CompositeMatrixTransform                                Self;
    typedef rstk::CachedMatrixTransform< TScalar, NDimensions >     Superclass;
    typedef itk::SmartPointer< Self >                               Pointer;
    typedef itk::SmartPointer< const Self >                         ConstPointer;

    itkTypeMacro( CompositeMatrixTransform, CachedMatrixTransform );
    itkNewMacro( Self );
    itkStaticConstMacro( Dimension, unsigned int, NDimensions );

    using typename Superclass::InterpolateModeType;

    typedef typename Superclass::ScalarType                          ScalarType;
    typedef typename Superclass::PointType                           PointType;
    typedef typename Superclass::VectorType                          VectorType;
    typedef typename Superclass::MatrixType                          MatrixType;

    /** Standard coordinate point type for this class. */
    typedef typename Superclass::InputPointType                      InputPointType;
    typedef typename Superclass::OutputPointType                     OutputPointType;

    /** Standard vector type for this class. */
    typedef typename Superclass::InputVectorType                     InputVectorType;
    typedef typename Superclass::OutputVectorType                    OutputVectorType;

    /** Standard covariant vector type for this class */
    typedef typename Superclass::InputCovariantVectorType            InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType           OutputCovariantVectorType;

    /** Standard vnl_vector type for this class. */
    typedef typename Superclass::InputVnlVectorType                  InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType                 OutputVnlVectorType;

    typedef typename Superclass::DisplacementFieldType               DisplacementFieldType;
    typedef typename DisplacementFieldType::Pointer                  DisplacementFieldPointer;
    typedef typename DisplacementFieldType::PixelType                DisplacementType;

    typedef rstk::SparseMatrixTransform< ScalarType >                TransformComponentType;
    typedef typename TransformComponentType::Pointer                 TransformComponentPointer;
    typedef typename TransformComponentType::ConstPointer            TransformComponentConstPointer;
    typedef std::vector<TransformComponentPointer>                   TransformsContainer;

    typedef typename Superclass::CoefficientsImageType               CoefficientsImageType;
    typedef typename Superclass::CoefficientsImageArray              CoefficientsImageArray;

    typedef typename Superclass::PointsList                          PointsList;
    typedef typename Superclass::DimensionVector                     DimensionVector;
    typedef typename Superclass::DimensionParameters                 DimensionParameters;
    typedef typename Superclass::DimensionParametersContainer        DimensionParametersContainer;

    typedef typename Superclass::FieldType                           FieldType;
    typedef typename FieldType::Pointer                              FieldPointer;
    typedef typename FieldType::ConstPointer                         FieldConstPointer;

    typedef rstk::BSplineSparseMatrixTransform< ScalarType >         BSplineComponentType;
    typedef typename BSplineComponentType::Pointer                   BSplineComponentPointer;

    typedef typename BSplineComponentType::AltCoeffType              AltCoeffType;
    typedef typename BSplineComponentType::AltCoeffPointer           AltCoeffPointer;


    itkSetMacro(NumberOfTransforms, size_t);
    itkGetConstMacro(NumberOfTransforms, size_t);

    void PushBackTransform(TransformComponentType* tf) {
    	this->m_Components.push_back(tf);
    	this->m_NumberOfTransforms = this->m_Components.size();
    }

    void PushBackCoefficients(const CoefficientsImageArray & images);
    void PushBackCoefficients(const FieldType* field);

    void Interpolate() override;
    void ComputeInverse() override {};
protected:
    CompositeMatrixTransform();
	~CompositeMatrixTransform(){};
    void PrintSelf( std::ostream& os, itk::Indent indent ) const override;

private:
	CompositeMatrixTransform( const Self & );
	void operator=( const Self & );

    void ComputeGrid();
    void ComputePoints();
    void ConsistencyCheck();

	TransformsContainer m_Components;
	size_t m_NumberOfTransforms;
};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "CompositeMatrixTransform.hxx"
#endif


#endif /* COMPOSITEMATRIXTRANSFORM_H_ */
