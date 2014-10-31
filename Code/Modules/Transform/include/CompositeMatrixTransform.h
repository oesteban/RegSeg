// --------------------------------------------------------------------------------------
// File:          CompositeMatrixTransform.h
// Date:          Oct 30, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
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

#ifndef COMPOSITEMATRIXTRANSFORM_H_
#define COMPOSITEMATRIXTRANSFORM_H_

#include <vector>
#include <iostream>
#include <itkDisplacementFieldTransform.h>
#include <itkMatrix.h>

#include "SparseMatrixTransform.h"
#include "BSplineSparseMatrixTransform.h"
#include "rstkMacro.h"

namespace rstk {

template< class TScalar, unsigned int NDimensions = 3u >
class CompositeMatrixTransform : public itk::DisplacementFieldTransform< TScalar, NDimensions >
{
public:
    /* Standard class typedefs. */
    typedef CompositeMatrixTransform          Self;
    typedef itk::DisplacementFieldTransform
    		         < TScalar, NDimensions > Superclass;
    typedef itk::SmartPointer< Self >         Pointer;
    typedef itk::SmartPointer< const Self >   ConstPointer;

    itkTypeMacro( CompositeMatrixTransform, DisplacementFieldTransform );
    itkNewMacro( Self );

    itkStaticConstMacro( Dimension, unsigned int, NDimensions );

    typedef typename Superclass::ScalarType          ScalarType;
    typedef itk::Point< ScalarType, Dimension >      PointType;
    typedef itk::Vector< ScalarType, Dimension >     VectorType;
    typedef itk::Matrix
    	    < ScalarType, Dimension, Dimension >     MatrixType;

    typedef itk::DisplacementFieldTransform< ScalarType, Dimension > DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer         DisplacementFieldTransformPointer;

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

    typedef BSplineSparseMatrixTransform< ScalarType >               TransformComponentType;
    typedef typename TransformComponentType::Pointer                 TransformComponentPointer;

    typedef typename TransformComponentType::FieldType               TransformFieldType;
    typedef typename TransformFieldType::Pointer                     TransformFieldPointer;
    typedef typename TransformFieldType::ConstPointer                TransformFieldConstPointer;
    typedef typename TransformComponentType::CoefficientsImageType   CoefficientsImageType;
    typedef typename TransformComponentType::CoefficientsImageArray  CoefficientsImageArray;
    typedef typename TransformComponentType::FieldType               CoefficientsFieldType;
    typedef typename CoefficientsFieldType::Pointer                  CoefficientsFieldTypePointer;
    typedef typename CoefficientsFieldType::ConstPointer             CoefficientsFieldTypeConstPointer;
    typedef typename TransformComponentType::DomainBase              DomainBaseType;
    typedef typename DomainBaseType::SizeType                        DomainSizeType;
    typedef typename DomainBaseType::DirectionType                   DomainDirectionType;
    typedef typename DomainBaseType::PointType                       DomainOriginType;
    typedef typename DomainBaseType::SpacingType                     DomainSpacingType;
    typedef typename std::vector<CoefficientsFieldTypeConstPointer>  CoefficientsContainer;

    itkSetMacro(NumberOfTransforms, size_t);
    itkGetConstMacro(NumberOfTransforms, size_t);

    void PushBackCoefficients(const CoefficientsFieldType* coeffs) {
    	this->m_Coefficients.push_back(coeffs);
    	this->m_NumberOfTransforms = this->m_Coefficients.size();
    }

    void SetPhysicalDomainInformation( const DomainBaseType* image );
    void Update() { this->Compute(); }

protected:
    CompositeMatrixTransform();
	~CompositeMatrixTransform(){};
    void PrintSelf( std::ostream& os, itk::Indent indent ) const;
    void Compute();

private:
	CompositeMatrixTransform( const Self & );
	void operator=( const Self & );

	CoefficientsContainer m_Coefficients;
	size_t m_NumberOfTransforms;
	DomainSizeType m_ReferenceSize;
	DomainDirectionType m_ReferenceDirection;
	DomainOriginType m_ReferenceOrigin;
	DomainSpacingType m_ReferenceSpacing;
};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "CompositeMatrixTransform.hxx"
#endif


#endif /* COMPOSITEMATRIXTRANSFORM_H_ */
