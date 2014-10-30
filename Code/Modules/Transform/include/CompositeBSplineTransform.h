// --------------------------------------------------------------------------------------
// File:          CompositeBSplineTransform.h
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

#ifndef COMPOSITEBSPLINETRANSFORM_H_
#define COMPOSITEBSPLINETRANSFORM_H_

#include <functional>

#include <itkTransform.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <itkMatrix.h>
#include <itkPointSet.h>
#include <itkImage.h>
#include <itkDefaultStaticMeshTraits.h>
#include <itkKernelFunctionBase.h>
#include <VNLSparseLUSolverTraits.h>
#include <itkDisplacementFieldTransform.h>
#include <itkImageHelper.h>
#include <itkImageTransformHelper.h>

#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_sparse_lu.h>

#include "rstkMacro.h"

namespace rstk {

template< class TScalar, unsigned int NDimensions = 3u >
class CompositeBSplineTransform: public itk::DisplacementFieldTransform< TScalar, NDimensions >
{
public:
    /* Standard class typedefs. */
    typedef CompositeBSplineTransform         Self;
    typedef itk::DisplacementFieldTransform< TScalar, NDimensions > Superclass;
    typedef itk::SmartPointer< Self >         Pointer;
    typedef itk::SmartPointer< const Self >   ConstPointer;

    itkTypeMacro( CompositeBSplineTransform, Transform );
    itkNewMacro( Self );

    itkStaticConstMacro( Dimension, unsigned int, NDimensions );

    typedef typename Superclass::ScalarType          ScalarType;
    typedef itk::Point< ScalarType, Dimension >      PointType;
    typedef itk::Vector< ScalarType, Dimension >     VectorType;
    typedef itk::Matrix
    	    < ScalarType, Dimension, Dimension >     MatrixType;

    typedef itk::DisplacementFieldTransform< ScalarType, Dimension > DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer         DisplacementFieldTransformPointer;

    typedef itk::FixedArray< DimensionVector, NDimensions > DimensionParametersContainer;

    typedef std::vector< PointType >                 PointsList;

    typedef itk::Matrix< ScalarType, Dimension, Dimension >        JacobianType;

    typedef itk::DefaultStaticMeshTraits<TScalar, NDimensions, NDimensions, TScalar, TScalar> PointSetTraitsType;
    typedef itk::PointSet<PointType, NDimensions, PointSetTraitsType>                                     PointSetType;
    typedef typename PointSetType::Pointer           PointSetPointer;


    /** Standard coordinate point type for this class. */
    typedef typename Superclass::InputPointType  InputPointType;
    typedef typename Superclass::OutputPointType OutputPointType;

    /** Standard vector type for this class. */
    typedef typename Superclass::InputVectorType  InputVectorType;
    typedef typename Superclass::OutputVectorType OutputVectorType;

    /** Standard covariant vector type for this class */
    typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;

    /** Standard vnl_vector type for this class. */
    typedef typename Superclass::InputVnlVectorType  InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType OutputVnlVectorType;

protected:
    CompositeBSplineTransform();
	~CompositeBSplineTransform(){};


private:
	CompositeBSplineTransform( const Self & );
	void operator=( const Self & );
};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "CompositeBSplineTransform.hxx"
#endif


#endif /* COMPOSITEBSPLINETRANSFORM_H_ */
