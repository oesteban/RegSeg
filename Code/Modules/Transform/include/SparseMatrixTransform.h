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
#include <itkImage.h>
#include <itkDefaultStaticMeshTraits.h>
#include <itkKernelFunctionBase.h>
#include <VNLSparseLUSolverTraits.h>
#include <itkDisplacementFieldTransform.h>

#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_sparse_lu.h>

#include "rstkMacro.h"

namespace rstk {

template< class TScalar, unsigned int NDimensions = 3u >
class SparseMatrixTransform: public itk::DisplacementFieldTransform< TScalar, NDimensions >
{
public:
	/* Standard class typedefs. */
	typedef SparseMatrixTransform             Self;
	typedef itk::DisplacementFieldTransform< TScalar, NDimensions > Superclass;
	typedef itk::SmartPointer< Self >         Pointer;
	typedef itk::SmartPointer< const Self >   ConstPointer;

	itkTypeMacro( SparseMatrixTransform, Transform );
	//itkCloneMacro(Self);
	itkNewMacro( Self );

	itkStaticConstMacro( Dimension, unsigned int, NDimensions );

	typedef typename Superclass::ScalarType          ScalarType;

	typedef itk::Point< ScalarType, Dimension >      PointType;
	typedef itk::Vector< ScalarType, Dimension >     VectorType;

    typedef itk::KernelFunctionBase<ScalarType>      KernelFunctionType;
    typedef typename KernelFunctionType::Pointer     KernelFunctionPointer;

    typedef VNLSparseLUSolverTraits< ScalarType >    SolverType;
	typedef typename SolverType::MatrixType          WeightsMatrix;
	typedef typename SolverType::VectorType          DimensionVector;
	typedef typename WeightsMatrix::row              SparseVectorType;

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

    typedef itk::FixedArray< ScalarType, itkGetStaticConstMacro(Dimension) >    ArrayType;

    typedef itk::Image< ScalarType, Dimension >                      CoefficientsImageType;
    typedef typename CoefficientsImageType::Pointer                  CoeffImagePointer;
    typedef typename CoefficientsImageType::ConstPointer             CoeffImageConstPointer;
    typedef itk::FixedArray< CoeffImagePointer, Dimension >          CoefficientsImageArray;

    /** Typedefs for specifying the extent of the grid. */
    typedef itk::ImageBase< Dimension >                   DomainBase;
    typedef typename DomainBase::Pointer                  DomainPointer;
    typedef itk::ImageRegion< Dimension >                 RegionType;
    typedef typename RegionType::IndexType                IndexType;
    typedef typename DomainBase::SizeType                 SizeType;
    typedef typename CoefficientsImageType::SpacingType   SpacingType;
    typedef typename CoefficientsImageType::DirectionType DirectionType;
    typedef typename CoefficientsImageType::PointType     OriginType;
    typedef typename CoefficientsImageType::SpacingType   PhysicalDimensionsType;
    typedef itk::ContinuousIndex< ScalarType, Dimension>  ContinuousIndexType;

    typedef itk::Image< VectorType, Dimension >           FieldType;
    typedef typename FieldType::Pointer                   FieldPointer;
    typedef typename FieldType::ConstPointer              FieldConstPointer;
    typedef typename std::vector< FieldPointer >          DerivativesType;

    /** Type of the input parameters. */
    typedef typename Superclass::ParametersType          ParametersType;
    typedef typename Superclass::ParametersValueType     ParametersValueType;

    /** Define the internal parameter helper used to access the field */
    typedef itk::ImageVectorOptimizerParametersHelper
    		                                 < ScalarType, Dimension, Dimension>
    OptimizerParametersHelperType;

	struct MatrixSectionType {
		WeightsMatrix *matrix;
		PointsList *vrows;
		PointsList *vcols;
		size_t section_id;
		size_t first_row;
		size_t num_rows;
	};

    itkSetMacro( ControlPointsSize, SizeType );
    itkGetConstMacro( ControlPointsSize, SizeType );

    itkSetObjectMacro(KernelFunction, KernelFunctionType);
    itkGetConstReferenceObjectMacro(KernelFunction, KernelFunctionType);

    itkGetConstMacro( NumberOfSamples, size_t );
    itkSetMacro( NumberOfSamples, size_t );

    itkGetConstMacro( NumberOfParameters, size_t );

    itkGetMacro( Derivatives, DerivativesType );

    void Interpolate();
    void UpdateField();
    //void ComputeCoeffDerivatives( void );
    void ComputeGradientField();
    void ComputeCoefficients();


	// Setters & Getters -----------------------------------------
	// Physical positions
	inline void SetOffGridPos  ( size_t id, const PointType pi );
	inline size_t AddOffGridPos  ( const PointType pi );
	void SetOffGridPositions( const PointsList points ) {
		this->m_OffGridPos = points;
		this->m_NumberOfSamples = points.size();
		this->Modified();
	}

	// Values off-grid (displacement vector of a node)
	inline bool       SetOffGridValue( const size_t id, VectorType pi );
	inline VectorType GetOffGridValue( const size_t id ) const;

	// Values on-grid (displacement vector of a grid point)
	//inline bool       SetOnGridValue ( const size_t id, VectorType pi );
	//inline VectorType GetOnGridValue ( const size_t id );

	// Coefficients of the interpolating kernels
	inline VectorType GetCoefficient ( const size_t id );

	virtual WeightsMatrix  GetPhi ();
	//itkGetConstMacro( Phi, WeightsMatrix );
	itkGetConstMacro( S, WeightsMatrix );
	itkGetConstMacro( OffGridValueMatrix, WeightsMatrix );

    void SetControlPointsSize( size_t s ) { this->m_ControlPointsSize.Fill( s ); }

    itkSetObjectMacro( Field, FieldType );
    itkGetConstObjectMacro( Field, FieldType );

    /** Get the array of coefficient images. */
    //const CoefficientsImageArray GetCoefficientsImages() const {
    //  return this->m_CoefficientsImages;
    //}
    itkGetConstMacro( CoefficientsImages, CoefficientsImageArray );
    void SetCoefficientsImages( const CoefficientsImageArray & images );
    void SetCoefficientsImage( size_t dim, const CoefficientsImageType* c );
    void SetCoefficientsVectorImage( const FieldType* f );

    void SetPhysicalDomainInformation( const DomainBase* image );
    void CopyGridInformation( const DomainBase* image );
    void SetOutputReference( const DomainBase* image );
    const FieldType* GetOutputField() const {
    	if (!this->m_UseImageOutput){
    		itkExceptionMacro( << "UseImageOutput is false, output field is not gridded" );
    	}
    	return this->m_OutputField;
    }

    /** Return the multithreader used by this class. */
    itk::MultiThreader * GetMultiThreader() const { return m_Threader; }

    itkSetClampMacro( NumberOfThreads, itk::ThreadIdType, 1, ITK_MAX_THREADS);
    itkGetConstReferenceMacro(NumberOfThreads, itk::ThreadIdType);
protected:
	SparseMatrixTransform();
	~SparseMatrixTransform(){};

	void ThreadedComputeMatrix( MatrixSectionType& section, itk::ThreadIdType threadId );
	itk::ThreadIdType SplitMatrixSection( itk::ThreadIdType i, itk::ThreadIdType num, MatrixSectionType& section );
	static ITK_THREAD_RETURN_TYPE ComputeThreaderCallback(void *arg);

	WeightsMatrix ComputeDerivativeMatrix( PointsList vrows, PointsList vcols, size_t dim );

	void InitializeCoefficientsImages();
	DimensionVector Vectorize( const CoefficientsImageType* image );
	WeightsMatrix VectorizeCoefficients();
	DimensionParametersContainer VectorizeField( const FieldType* image );
	WeightsMatrix MatrixField( const FieldType* image );

	inline ScalarType Evaluate( const VectorType r ) const;
	inline ScalarType EvaluateDerivative( const VectorType r, size_t dim ) const;
	//virtual size_t GetSupport() const = 0;

	/* Field domain definitions */
	SizeType               m_ControlPointsSize;
	SpacingType            m_ControlPointsSpacing;
	OriginType             m_ControlPointsOrigin;
	DirectionType          m_ControlPointsDirection;
	DirectionType          m_ControlPointsDirectionInverse;

	size_t                 m_NumberOfSamples;     // This is N mesh points
	size_t                 m_NumberOfParameters;  // This is K nodes

	PointsList m_OffGridPos;           // m_N points in the mesh
	PointsList m_OnGridPos;

	WeightsMatrix m_OffGridValueMatrix;     // m_N points in the mesh
	CoefficientsImageArray m_CoefficientsImages;


	//DimensionParametersContainer m_CoeffDerivative;  // Serialized k values in a grid
	//DimensionVector m_Jacobian[Dimension][Dimension]; // Serialized k dimxdim matrices in a grid

	WeightsMatrix   m_Phi;
	WeightsMatrix   m_S;
	WeightsMatrix   m_SPrime[Dimension];


	bool            m_GridDataChanged;
	bool            m_ControlDataChanged;
	bool            m_UseImageOutput;

	KernelFunctionPointer m_KernelFunction;
	KernelFunctionPointer m_DerivativeKernel;
	KernelFunctionPointer m_SecondDerivativeKernel;
	FieldPointer          m_Field;
	DerivativesType       m_Derivatives;
	FieldPointer          m_OutputField;

private:
	SparseMatrixTransform( const Self & );
	void operator=( const Self & );

	enum MatrixType { PHI, S };

	struct SMTStruct {
		SparseMatrixTransform *Transform;
		MatrixType type;
		WeightsMatrix *matrix;
		PointsList *vrows;
		PointsList *vcols;
	};

	void ComputeMatrix( MatrixType type );

	/** Support processing data in multiple threads. */
	itk::MultiThreader::Pointer m_Threader;
	itk::ThreadIdType           m_NumberOfThreads;
};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SparseMatrixTransform.hxx"
#endif
#endif /* SPARSEMATRIXTRANSFORM_H_ */
