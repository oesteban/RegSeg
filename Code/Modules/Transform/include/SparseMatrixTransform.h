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
#include <itkMatrix.h>
#include <itkPointSet.h>
#include <itkDefaultStaticMeshTraits.h>
#include <itkKernelFunctionBase.h>

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
	//itkNewMacro( Self );

	itkStaticConstMacro( Dimension, unsigned int, NDimensions );

	typedef typename Superclass::ScalarType ScalarType;

	typedef itk::Point< ScalarType, Dimension >      PointType;
	typedef itk::Vector< ScalarType, Dimension >     VectorType;

    typedef itk::KernelFunctionBase<ScalarType>      KernelFunctionType;

	typedef vnl_sparse_matrix< ScalarType >          WeightsMatrix;
	typedef vnl_vector< ScalarType >                 DimensionVector;
	typedef itk::FixedArray< DimensionVector, NDimensions > ParametersType;

	typedef std::vector< PointType >                 PointsList;

	typedef itk::Matrix< ScalarType, Dimension, Dimension >        JacobianType;

	typedef itk::DefaultStaticMeshTraits<TScalarType, NDimensions, NDimensions, TScalarType, TScalarType> PointSetTraitsType;
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

    typedef itk::FixedArray< ScalarType, itkGetStaticConstMacro(Dimension) >    ArrayType;

    itkSetObjectMacro(KernelFunction, KernelFunctionType);
    itkGetConstReferenceObjectMacro(KernelFunction, KernelFunctionType);

    void SetN( size_t N );
    itkGetConstMacro(N, size_t );

    void SetNumberOfParameters( size_t K );
    itkGetConstMacro(NumberOfParameters, size_t );

    void Interpolate( void );
    void ComputeOnGridValues( void );
    void ComputeCoeffDerivatives( void );
    void ComputeJacobian( void );
    void ComputeCoefficients();


	// Setters & Getters -----------------------------------------
	// Physical positions
	inline void SetOffGridPos  ( size_t id, const PointType pi );
	inline void SetOnGridPos   ( size_t id, const PointType pi );

	// Values off-grid (displacement vector of a node)
	inline bool       SetOffGridValue( const size_t id, VectorType pi );
	inline VectorType GetOffGridValue( const size_t id ) const;

	// Values on-grid (displacement vector of a grid point)
	inline bool       SetOnGridValue ( const size_t id, VectorType pi );
	inline VectorType GetOnGridValue ( const size_t id );

	// Coefficients of the interpolating kernels
	inline bool       SetCoefficient ( const size_t id, VectorType pi );
	inline VectorType GetCoefficient ( const size_t id );

	// Derivative of coefficients
	inline VectorType GetCoeffDerivative ( const size_t id );

	// Sparse matrix of kernel weights
	//inline VectorType GetOffGridWeight ( const size_t id );
	//inline VectorType GetOnGridWeight  ( const size_t id );

	inline JacobianType GetJacobian    ( const size_t id );

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

    const WeightsMatrix * GetPhi() const { return this->m_Phi; }

    typedef itk::Image< ScalarType, Dimension >                      CoefficientsImageType;
    typedef typename CoefficientsImageType::Pointer                  CoeffImagePointer;
    typedef itk::FixedArray< CoeffImagePointer, Dimension >          CoefficientImageArray;

    void SetCoefficientImages( const CoefficientImageArray & images );
    /** Get the array of coefficient images. */
    const CoefficientImageArray GetCoefficientImages() const {
      return this->m_CoefficientImages;
    }

    /** Typedefs for specifying the extent of the grid. */
    typedef itk::ImageBase< Dimension >                   DomainBase;
    typedef typename DomainBase::Pointer                  DomainPointer;
    typedef itk::ImageRegion< Dimension >                 RegionType;
    typedef typename RegionType::IndexType                IndexType;
    typedef typename RegionType::SizeType                 SizeType;
    typedef typename CoefficientsImageType::SpacingType   SpacingType;
    typedef typename CoefficientsImageType::DirectionType DirectionType;
    typedef typename CoefficientsImageType::PointType     OriginType;
    typedef typename CoefficientsImageType::SpacingType   PhysicalDimensionsType;
    typedef itk::ContinuousIndex< ScalarType, Dimension>  ContinuousIndexType;

    itkSetMacro( ControlPointsSize, SizeType );
    itkGetConstMacro( ControlPointsSize, SizeType );

    void SetPhysicalDomainInformation( const DomainBase* image );
protected:
	SparseMatrixTransform();
	~SparseMatrixTransform(){};


	void ComputePhi();
	void ComputeS();
	void ComputeSPrime();

	virtual ScalarType Evaluate( VectorType r ) = 0;

	/* Field domain definitions */
	SizeType               m_ControlPointsSize;
	SpacingType            m_ControlPointsSpacing;
	OriginType             m_ControlPointsOrigin;
	DirectionType          m_ControlPointsDirection;
	DirectionType          m_ControlPointsDirectionInverse;

	size_t                 m_NumberOfSamples;     // This is N mesh points
	size_t                 m_NumberOfParameters;  // This is K nodes


	PointsList m_OffGridPos;           // m_N points in the mesh

	ParametersType m_OffGridValue;     // m_N points in the mesh
	ParametersType m_OnGridValue;      // Serialized k values in a grid
	ParametersType m_Coeff;            // Serialized k coefficients in the grid
	ParametersType m_CoeffDerivative;  // Serialized k values in a grid
	CoefficientImageArray m_CoefficientImages;

	DimensionVector m_Jacobian[Dimension][Dimension]; // Serialized k dimxdim matrices in a grid

	WeightsMatrix   m_Phi;
	WeightsMatrix   m_S;
	WeightsMatrix   m_SPrime[Dimension];


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
