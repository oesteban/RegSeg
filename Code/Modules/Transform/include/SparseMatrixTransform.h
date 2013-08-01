// --------------------------------------------------------------------------------------
// File:          SparseMatrixTransform.h
// Date:          Jun 6, 2013
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

#ifndef SPARSEMATRIXTRANSFORM_H_
#define SPARSEMATRIXTRANSFORM_H_

#include <itkTransform.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <itkPointSet.h>
#include <itkDefaultStaticMeshTraits.h>
#include <itkKernelFunctionBase.h>

#include "RadialBasisFunction.h"
#include "VectorIDWBasisFunction.h"

#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_sparse_lu.h>

namespace rstk {

template< class TScalarType = double, unsigned int NDimensions = 3u >
class SparseMatrixTransform: public itk::Transform< TScalarType, NDimensions, NDimensions >
{
public:
	/* Standard class typedefs. */
	typedef SparseMatrixTransform             Self;
	typedef itk::Transform< TScalarType, NDimensions, NDimensions > Superclass;
	typedef itk::SmartPointer< Self >         Pointer;
	typedef itk::SmartPointer< const Self >   ConstPointer;

	itkTypeMacro( SparseMatrixTransform, Transform );
	itkNewMacro( Self );

	itkStaticConstMacro( Dimension, unsigned int, NDimensions );

	typedef typename Superclass::ScalarType ScalarType;
	typedef typename Superclass::ParametersType ParametersType;

	typedef itk::Point< ScalarType, Dimension >      PointType;
	typedef itk::Vector< ScalarType, Dimension >     VectorType;

    typedef itk::KernelFunctionBase<ScalarType>      KernelFunctionType;


	typedef RBF::RadialBasisFunction
			< PointType, TScalarType, NDimensions >  RBFType;

	typedef RBF::VectorIDWBasisFunction
			< PointType, TScalarType, NDimensions >  DefaultRBFType;

	typedef vnl_sparse_matrix< ScalarType >          WeightsMatrix;
	typedef vnl_vector< ScalarType >                 DimensionVector;

	typedef std::vector< PointType >                 PointsList;


	typedef itk::DefaultStaticMeshTraits<TScalarType, NDimensions, NDimensions, TScalarType, TScalarType> PointSetTraitsType;
	typedef itk::PointSet<PointType, NDimensions, PointSetTraitsType>                                     PointSetType;
	typedef typename PointSetType::Pointer           PointSetPointer;

	typedef typename Superclass::JacobianType JacobianType;

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

    typedef itk::FixedArray< ScalarType, itkGetStaticConstMacro(Dimension) >    ArrayType;


    itkSetObjectMacro(KernelFunction, KernelFunctionType);
    itkGetConstReferenceObjectMacro(KernelFunction, KernelFunctionType);

    itkSetMacro(Sigma, ArrayType);
    itkGetConstReferenceMacro(Sigma, ArrayType);

    void SetN( size_t N );
    void SetK( size_t K );
    void SetNumberOfParameters( size_t N, size_t K ) {
    	this->SetN( N );
    	this->SetK( K );
    }

    void ComputeWeights( void );
    void ComputeNodesData( void );
	void Interpolate( void );

	inline void SetPoint( size_t id, const PointType pi );
	inline void SetNode ( size_t id, const PointType pi );

	inline VectorType GetPointData  ( const size_t id );
	inline VectorType GetNodeData   ( const size_t id );
	inline VectorType GetNodeWeight ( const size_t id );
	inline VectorType GetCoefficient( const size_t id );

	inline bool SetPointData( const size_t id, VectorType pi );
	inline bool SetNodeData( const size_t id, VectorType pi );
	inline bool SetCoefficient( const size_t id, VectorType pi );

	// Virtual members inherited from Transform
	virtual void SetParameters(const ParametersType & parameters);

	virtual void SetFixedParameters(const ParametersType &)                       \
		    {                                                                                             \
		      itkExceptionMacro(                                                                          \
		        << "TransformVector(const InputVectorType &) is not implemented for KernelTransform");    \
		    }

	virtual OutputPointType TransformPoint(const InputPointType  & point) const
	{
	  return point;
	}

	/** These vector transforms are not implemented for this transform */
    using Superclass::TransformVector;
    virtual OutputVectorType TransformVector(const InputVectorType &) const                       \
    {                                                                                             \
      itkExceptionMacro(                                                                          \
        << "TransformVector(const InputVectorType &) is not implemented for KernelTransform");    \
    }

    virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const                 \
    {                                                                                             \
      itkExceptionMacro(                                                                          \
        << "TransformVector(const InputVnlVectorType &) is not implemented for KernelTransform"); \
    }

    /**  Method to transform a CovariantVector. */
    using Superclass::TransformCovariantVector;
    virtual OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType &) const           \
    {                                                                                                            \
      itkExceptionMacro(                                                                                         \
        << "TransformCovariantVector(const InputCovariantVectorType &) is not implemented for KernelTransform"); \
    }
    /** Compute the Jacobian Matrix of the transformation at one point */
    virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const           \
    {                                                                                 \
      itkExceptionMacro( "ComputeJacobianWithRespectToPosition not yet implemented "  \
                         "for " << this->GetNameOfClass() );                          \
    }

    virtual void ComputeJacobianWithRespectToPosition(const InputPointType &,
                                                      JacobianType &) const           \
    {                                                                                 \
      itkExceptionMacro( "ComputeJacobianWithRespectToPosition not yet implemented "  \
                         "for " << this->GetNameOfClass() );                          \
    }


//    void SetRadialBasisFunction( RBFType* rbf ) {
//    	this->m_RadialBasisFunction = rbf;
//    }

protected:
	SparseMatrixTransform();
	~SparseMatrixTransform(){};


	void ComputePhi( void );
	void ComputeS( void );

	PointsList m_Points; // Nc points in the mesh
	PointsList m_Nodes;    // Serialized k points in a grid

	DimensionVector m_PointsData[3];     // Nc points in the mesh
	DimensionVector m_NodesData[3];      // Serialized k values in a grid
	DimensionVector m_Coeff[3];          // Serialized k coefficients in the grid
	DimensionVector m_TempNodesData[3];  // Serialized k values in a grid

	WeightsMatrix   m_Phi;
	WeightsMatrix   m_S;
	WeightsMatrix   m_InvertPhi;
	size_t          m_N;
	size_t          m_K;

	vnl_sparse_lu* m_System;

	typename KernelFunctionType::Pointer  m_KernelFunction;
	ScalarType m_KernelNorm;
	ArrayType m_Sigma;

	bool            m_GridDataChanged;
	bool            m_ControlDataChanged;
private:
	SparseMatrixTransform( const Self & );
	void operator=( const Self & );

};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SparseMatrixTransform.hxx"
#endif
#endif /* SPARSEMATRIXTRANSFORM_H_ */
