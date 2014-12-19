// --------------------------------------------------------------------------------------
// File:          CachedMatrixTransform.h
// Date:          Nov 25, 2014
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

#ifndef CACHEDMATRIXTRANSFORM_H_
#define CACHEDMATRIXTRANSFORM_H_

#include <functional>

#include <itkTransform.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <itkMatrix.h>
#include <itkImage.h>
#include <itkDefaultStaticMeshTraits.h>
#include <itkKernelFunctionBase.h>
#include <VNLSparseLUSolverTraits.h>
#include <itkImageHelper.h>
#include <itkImageTransformHelper.h>
#include <itkInvertDisplacementFieldImageFilter.h>

#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_sparse_lu.h>

#include "RBFFieldTransform.h"
#include "rstkMacro.h"

namespace rstk {

template< class TScalar, unsigned int NDimensions = 3u >
class CachedMatrixTransform: public RBFFieldTransform< TScalar, NDimensions >
{
public:
    /* Standard class typedefs. */
    typedef CachedMatrixTransform                     Self;
    typedef RBFFieldTransform< TScalar, NDimensions > Superclass;
    typedef itk::SmartPointer< Self >                 Pointer;
    typedef itk::SmartPointer< const Self >           ConstPointer;
    
    itkTypeMacro( CachedMatrixTransform, RBFFieldTransform );
    itkStaticConstMacro( Dimension, unsigned int, NDimensions );
    
    typedef enum {
      UNKNOWN,
      GRID_MODE,
      POINTS_MODE
    } InterpolateModeType;

    
    typedef typename Superclass::ScalarType                                     ScalarType;
    typedef itk::Point< ScalarType, Dimension >                                 PointType;
    typedef itk::Vector< ScalarType, Dimension >                                VectorType;
    typedef itk::Matrix
    	    < ScalarType, Dimension, Dimension >                                MatrixType;

    typedef vnl_sparse_matrix< ScalarType >                                     WeightsMatrix;
    typedef typename WeightsMatrix::row                                         SparseMatrixRowType;
    typedef vnl_vector< ScalarType >                                            DimensionVector;
    typedef vnl_matrix< ScalarType >                                            DimensionMatrixType;

    typedef VNLSparseLUSolverTraits< double >                                   SolverTypeTraits;
    typedef typename SolverTypeTraits::SolverType                               SolverType;
    typedef typename SolverTypeTraits::MatrixType                               SolverMatrix;
    typedef typename SolverTypeTraits::VectorType                               SolverVector;
    typedef typename SolverMatrix::pair_t                                       SolverPair;

    typedef itk::FixedArray< DimensionVector, NDimensions >                     DimensionParameters;
    typedef itk::FixedArray< DimensionVector*, NDimensions >                    DimensionParametersContainer;

    typedef std::vector< PointType >                                            PointsList;

    typedef itk::Matrix< ScalarType, Dimension, Dimension >                     JacobianType;
    
    // typedef itk::DefaultStaticMeshTraits
    // 		              <TScalar, NDimensions, NDimensions, TScalar, TScalar> PointSetTraitsType;
    // typedef itk::PointSet<PointType, NDimensions, PointSetTraitsType>           PointSetType;
    // typedef typename PointSetType::Pointer                                      PointSetPointer;


    /** Standard coordinate point type for this class. */
    typedef typename Superclass::InputPointType                                 InputPointType;
    typedef typename Superclass::OutputPointType                                OutputPointType;

    /** Standard vector type for this class. */
    typedef typename Superclass::InputVectorType                                InputVectorType;
    typedef typename Superclass::OutputVectorType                               OutputVectorType;

    /** Standard covariant vector type for this class */
    typedef typename Superclass::InputCovariantVectorType                       InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType                      OutputCovariantVectorType;

    /** Standard vnl_vector type for this class. */
    typedef typename Superclass::InputVnlVectorType                             InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType                            OutputVnlVectorType;

    typedef itk::FixedArray< ScalarType, itkGetStaticConstMacro(Dimension) >    ArrayType;
    typedef itk::FixedArray< PointType, 2 >                                     ExtentType;

    typedef itk::Image< ScalarType, Dimension >                                 CoefficientsImageType;
    typedef typename CoefficientsImageType::Pointer                             CoeffImagePointer;
    typedef typename CoefficientsImageType::ConstPointer                        CoeffImageConstPointer;
    typedef itk::FixedArray< CoeffImagePointer, Dimension >                     CoefficientsImageArray;

    /** Typedefs for specifying the extent of the grid. */
    typedef itk::ImageBase< Dimension >                                         DomainBase;
    typedef typename DomainBase::Pointer                                        DomainPointer;
    typedef itk::ImageRegion< Dimension >                                       RegionType;
    typedef typename RegionType::IndexType                                      IndexType;
    typedef typename DomainBase::SizeType                                       SizeType;
    typedef typename CoefficientsImageType::SpacingType                         SpacingType;
    typedef typename CoefficientsImageType::DirectionType                       DirectionType;
    typedef typename CoefficientsImageType::PointType                           OriginType;
    typedef typename CoefficientsImageType::SpacingType                         PhysicalDimensionsType;
    typedef itk::ContinuousIndex< ScalarType, Dimension>                        ContinuousIndexType;
    typedef typename IndexType::OffsetType                                      OffsetType;
    typedef typename OffsetType::OffsetValueType                                OffsetValueType;
    typedef OffsetValueType                                                     OffsetTableType[Dimension+1];

    typedef itk::Image< VectorType, Dimension >                                 FieldType;
    typedef typename FieldType::Pointer                                         FieldPointer;
    typedef typename FieldType::ConstPointer                                    FieldConstPointer;
    //typedef typename std::vector< FieldPointer >                                DerivativesType;


    typedef itk::InvertDisplacementFieldImageFilter<FieldType, FieldType>       InvertFieldFilter;
    typedef typename InvertFieldFilter::Pointer                                 InvertFieldPointer;

    /** Type of the input parameters. */
    typedef typename Superclass::ParametersType                                 ParametersType;
    typedef typename Superclass::ParametersValueType                            ParametersValueType;
    typedef typename Superclass::LightObject                                    LightObject;

    /** Define the internal parameter helper used to access the field */
    typedef itk::ImageVectorOptimizerParametersHelper
    		      < ScalarType, Dimension, Dimension>    OptimizerParametersHelperType;
    typedef itk::ImageHelper< Dimension, Dimension >     Helper;
    typedef itk::ImageTransformHelper< Dimension, Dimension - 1, Dimension - 1, ScalarType, ScalarType > TransformHelper;

    typedef size_t                                                              PointIdentifier;
    typedef std::vector< PointIdentifier >                                      PointIdContainer;


    itkGetConstMacro( InterpolationMode, InterpolateModeType );
    itkGetConstMacro( NumberOfPoints, size_t );

    // Physical domain information
    itkSetMacro( DomainExtent, ExtentType );
    itkGetConstMacro( DomainExtent, ExtentType );
    void SetDomainExtent( const DomainBase* image );

    itkGetConstMacro( PointLocations, PointsList );
	itkGetConstMacro( PointValues, DimensionParameters );
    inline VectorType GetPointValue( const size_t id ) const;


	// Physical positions, will define interpolation mode
    void SetOutputReference( const DomainBase* image );
	void SetOutputPoints( const PointsList points );
	void SetOutputPoints( const PointsList points, const PointIdContainer valid);

    virtual void Interpolate() = 0;
    virtual void ComputeInverse() = 0;

    virtual void SetFixedParametersFromImage(const DomainBase* image);
protected:
	CachedMatrixTransform();
	~CachedMatrixTransform(){};

	DimensionVector Vectorize( const CoefficientsImageType* image );
	DimensionParameters VectorizeField( const FieldType* image );

	// Setters & Getters -----------------------------------------
	inline void SetPointLocation ( size_t id, const PointType pi ) {
		if ( id >= this->m_NumberOfPoints ) {
			itkExceptionMacro(<< "Trying to set sample with id " << id << ", when NumberOfPoints is set to " << this->m_NumberOfPoints );
		}
		if ( this->m_PointLocations.size() == 0 ) {
			this->m_PointLocations.resize( this->m_NumberOfPoints );
		}
		this->m_PointLocations[id] = pi;
	}

	inline size_t AddPointLocation  ( const PointType pi ) {
		this->m_PointLocations.push_back( pi );
		this->m_NumberOfPoints++;
		return (this->m_NumberOfPoints-1);
	}


	ExtentType 					 m_DomainExtent;
	DirectionType                m_ReferenceDirection;
	SpacingType                  m_ReferenceSpacing;
	SizeType                     m_ReferenceSize;
	PointType                    m_ReferenceOrigin;


	PointsList                   m_PointLocations;     // m_N points in the mesh
	PointIdContainer             m_ValidLocations;
	size_t                       m_NumberOfPoints;     // This is N mesh points
	DimensionParameters          m_PointValues;     // m_N points in the mesh

	bool                         m_UseImageOutput;
	InterpolateModeType          m_InterpolationMode;
private:
	CachedMatrixTransform( const Self & );
	void operator=( const Self & );
};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "CachedMatrixTransform.hxx"
#endif
#endif /* CACHEDMATRIXTRANSFORM_H_ */
