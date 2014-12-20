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

#include <functional>

#include "CachedMatrixTransform.h"
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
class SparseMatrixTransform: public rstk::CachedMatrixTransform< TScalar, NDimensions >
{
public:
    /* Standard class typedefs. */
    typedef SparseMatrixTransform                                  Self;
    typedef rstk::CachedMatrixTransform< TScalar, NDimensions >    Superclass;
    typedef itk::SmartPointer< Self >                              Pointer;
    typedef itk::SmartPointer< const Self >                        ConstPointer;
    
    itkStaticConstMacro( Dimension, unsigned int, NDimensions );
    itkTypeMacro( SparseMatrixTransform, CachedMatrixTransform );
    itkNewMacro( Self );
    
    using Superclass::InterpolateModeType;

    typedef typename Superclass::ScalarType                        ScalarType;
    typedef typename Superclass::PointType                         PointType;
    typedef typename Superclass::VectorType                        VectorType;
    typedef typename Superclass::MatrixType                        MatrixType;

    typedef typename Superclass::WeightsMatrix                     WeightsMatrix;
    typedef typename Superclass::SparseMatrixRowType               SparseMatrixRowType;
    typedef typename Superclass::DimensionVector                   DimensionVector;
    typedef typename Superclass::DimensionMatrixType               DimensionMatrixType;

    typedef typename Superclass::SolverTypeTraits                  SolverTypeTraits;
    typedef typename Superclass::SolverType                        SolverType;
    typedef typename Superclass::SolverMatrix                      SolverMatrix;
    typedef typename Superclass::SolverVector                      SolverVector;
    typedef typename Superclass::SolverPair                        SolverPair;
    
    typedef typename Superclass::DimensionParameters               DimensionParameters;
    typedef typename Superclass::DimensionParametersContainer      DimensionParametersContainer;
    typedef typename Superclass::PointsList                        PointsList;
    typedef typename Superclass::JacobianType                      JacobianType;
    // typedef typename Superclass::PointSetTraitsType                PointSetTraitsType;
    // typedef typename Superclass::PointSetType                      PointSetType;
    // typedef typename Superclass::PointSetPointer                   PointSetPointer;


    /** Standard coordinate point type for this class. */
    typedef typename Superclass::InputPointType                    InputPointType;
    typedef typename Superclass::OutputPointType                   OutputPointType;

    /** Standard vector type for this class. */
    typedef typename Superclass::InputVectorType                   InputVectorType;
    typedef typename Superclass::OutputVectorType                  OutputVectorType;

    /** Standard covariant vector type for this class */
    typedef typename Superclass::InputCovariantVectorType          InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType         OutputCovariantVectorType;

    /** Standard vnl_vector type for this class. */
    typedef typename Superclass::InputVnlVectorType                InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType               OutputVnlVectorType;

    typedef typename Superclass::ArrayType                         ArrayType;
    typedef typename Superclass::CoefficientsImageType             CoefficientsImageType;
    typedef typename Superclass::CoeffImagePointer                 CoeffImagePointer;
    typedef typename Superclass::CoeffImageConstPointer            CoeffImageConstPointer;
    typedef typename Superclass::CoefficientsImageArray            CoefficientsImageArray;

    /** Typedefs for specifying the extent of the grid. */
    typedef typename Superclass::DomainBase                        DomainBase;
    typedef typename Superclass::DomainPointer                     DomainPointer;
    typedef typename Superclass::RegionType                        RegionType;
    typedef typename Superclass::IndexType                         IndexType;
    typedef typename Superclass::SizeType                          SizeType;
    typedef typename Superclass::SpacingType                       SpacingType;
    typedef typename Superclass::DirectionType                     DirectionType;
    typedef typename Superclass::OriginType                        OriginType;
    typedef typename Superclass::PhysicalDimensionsType            PhysicalDimensionsType;
    typedef typename Superclass::ContinuousIndexType               ContinuousIndexType;
    typedef typename Superclass::OffsetType                        OffsetType;
    typedef typename Superclass::OffsetValueType                   OffsetValueType;
    typedef typename Superclass::OffsetTableType                   OffsetTableType;

    typedef typename Superclass::FieldType                         FieldType;
    typedef typename Superclass::FieldPointer                      FieldPointer;
    typedef typename Superclass::FieldConstPointer                 FieldConstPointer;
    typedef typename Superclass::InvertFieldFilter                 InvertFieldFilter;
    typedef typename Superclass::InvertFieldPointer                InvertFieldPointer;
    //typedef typename std::vector< FieldPointer >                     DerivativesType;

    /** Type of the input parameters. */
    typedef typename Superclass::ParametersType            ParametersType;
    typedef typename Superclass::ParametersValueType       ParametersValueType;

    typedef itk::PointSet<OutputVectorType, Dimension>                          AltCoeffType;
    typedef typename AltCoeffType::Pointer                                      AltCoeffPointer;
    typedef typename AltCoeffType::PointsContainerPointer                       AltCoeffContainerPointer;
    typedef typename AltCoeffType::PointDataContainerPointer                    AltCoeffDataPointer;

    /** Define the internal parameter helper used to access the field */
    typedef typename Superclass::OptimizerParametersHelperType OptimizerParametersHelperType;
    typedef typename Superclass::Helper                        Helper;
    typedef typename Superclass::TransformHelper               TransformHelper;

    typedef typename Superclass::PointIdContainer              PointIdContainer;


    typedef itk::KernelFunctionBase<ScalarType>      KernelFunctionType;
    typedef typename KernelFunctionType::Pointer     KernelFunctionPointer;

    itkSetObjectMacro(KernelFunction, KernelFunctionType);
    itkGetConstReferenceObjectMacro(KernelFunction, KernelFunctionType);

    itkGetConstMacro(ControlGridSize, SizeType);
    itkSetMacro(ControlGridSize, SizeType);
    itkGetConstMacro(ControlGridSpacing, SpacingType );
    itkSetMacro(ControlGridSpacing, SpacingType );
    itkGetConstMacro(ControlGridOrigin, PointType );
    itkSetMacro(ControlGridOrigin, PointType );

    itkGetConstMacro(MaximumDisplacement, SpacingType);

    void SetControlGridSize( size_t s ) {
    	SizeType size; size.Fill(s);
    	this->SetControlGridSize(size);
    }

    void SetControlGridSpacing( float s ) {
    	SpacingType ss;
        ss.Fill(s);
    	this->SetControlGridSpacing(ss);
    }

    void SetControlGridSpacing( itk::FixedArray< double, Dimension > s) {
    	SpacingType ss;
        for( size_t i = 0; i<Dimension; i++ )
            ss[i] = s[i];
    	this->SetControlGridSpacing(ss);
    }

    void SetNumberOfPoints(size_t n) {
    	if( this->m_InterpolationMode == Superclass::GRID_MODE ) {
    		itkExceptionMacro(<< "SetNumberOfSamples should not be used with OutputReference");
    	}

    	size_t oldsize = this->m_PointValues[0].size();
    	if ( oldsize!=0 && oldsize!=n ) {
    		itkWarningMacro( << "this action will empty off-grid values set so-far" );
    	}

    	this->m_NumberOfPoints = n;

    	for (size_t i=0; i<Dimension; i++) {
    		this->m_PointValues[i].set_size( n );
    		this->m_PointValues[i].fill( 0.0 );
    	}

    	this->Modified();
    }

    itkGetMacro( Derivatives, CoefficientsImageArray );
    itkGetConstObjectMacro( GradientField, FieldType );

    void InterpolateGradient() { this->Interpolate( this->VectorizeDerivatives() ); };
    void UpdateField() { this->UpdateField( this->VectorizeCoefficients() ); }
    void ComputeInverse();

    //void ComputeCoeffDerivatives( void );
    void ComputeGradientField();
    void ComputeCoefficients();

	// Values off-grid (displacement vector of a node)
	inline bool       SetPointValue( const size_t id, VectorType pi );

	// Coefficients of the interpolating kernels
	inline VectorType GetCoefficient ( const size_t id );

	virtual const WeightsMatrix*  GetPhi (const bool onlyvalid = true);
	virtual const WeightsMatrix*  GetS (){
		return &this->m_S;
	}

    itkSetObjectMacro( Field, FieldType );
    itkGetConstObjectMacro( Field, FieldType );

    /** Get the array of coefficient images. */
    itkGetConstMacro( CoefficientsImages, CoefficientsImageArray );

    void SetCoefficientsImages( const CoefficientsImageArray & images );
    void SetCoefficientsImage( size_t dim, const CoefficientsImageType* c );
    void SetCoefficientsVectorImage( const FieldType* f );
    void AddCoefficientsVectorImage( const FieldType* f );
    const FieldType* GetCoefficientsVectorImage();


    void SetControlGridInformation( const DomainBase* image );

	void Initialize();
	void Interpolate() { this->Interpolate( this->VectorizeCoefficients() ); }
	AltCoeffPointer GetFlatParameters();

    /** Return the multithreader used by this class. */
    itk::MultiThreader * GetMultiThreader() const { return m_Threader; }
    itkSetClampMacro( NumberOfThreads, itk::ThreadIdType, 1, ITK_MAX_THREADS);
    itkGetConstReferenceMacro(NumberOfThreads, itk::ThreadIdType);
protected:
	SparseMatrixTransform();
	~SparseMatrixTransform(){};

	enum WeightsMatrixType { PHI, S, SPRIME, PHI_INV };

	struct MatrixSectionType {
		WeightsMatrix *matrix;
		PointsList *vrows;
		PointsList *vcols;
		size_t section_id;
		size_t first_row;
		size_t num_rows;
		size_t dim;
	};

	typedef ScalarType (Self::*FunctionalCallback)( const VectorType, const size_t );

	struct SMTStruct {
		SparseMatrixTransform *Transform;
		WeightsMatrixType type;
		WeightsMatrix* matrix;
		size_t dim;
		PointsList *vrows;
		PointsList *vcols;
	};

	void Interpolate( const DimensionParameters& coeff );
	void UpdateField( const DimensionParameters& coeff );
	void InvertPhi();

	void ThreadedComputeMatrix( MatrixSectionType& section, FunctionalCallback func, itk::ThreadIdType threadId );
	itk::ThreadIdType SplitMatrixSection( itk::ThreadIdType i, itk::ThreadIdType num, MatrixSectionType& section );
	static ITK_THREAD_RETURN_TYPE ComputeThreaderCallback(void *arg);

	void InitializeCoefficientsImages();
	DimensionVector Vectorize( const CoefficientsImageType* image );
	//WeightsMatrix VectorizeCoefficients();
	DimensionParameters VectorizeCoefficients() const;
	DimensionParameters VectorizeDerivatives() const;
	DimensionParameters VectorizeField( const FieldType* image );
	WeightsMatrix MatrixField( const FieldType* image );

	inline ScalarType EvaluateKernel( const VectorType r, const size_t dim = 0 );
	inline ScalarType EvaluateDerivative( const VectorType r, const size_t dim );


	/* Field domain definitions */
	SizeType                     m_ControlGridSize;
	SpacingType                  m_ControlGridSpacing;
	PointType                    m_ControlGridOrigin;
	DirectionType                m_ControlGridDirection;
	DirectionType                m_ControlGridDirectionInverse;
	MatrixType                   m_ControlGridIndexToPhysicalPoint;
	MatrixType                   m_ControlGridPhysicalPointToIndex;
	SpacingType                  m_MaximumDisplacement;

	AltCoeffPointer              m_FlatCoeffs;

	CoefficientsImageArray       m_CoefficientsImages;
	CoefficientsImageArray       m_Derivatives;

	//DimensionParameters m_CoeffDerivative;  // Serialized k values in a grid
	//DimensionVector m_Jacobian[Dimension][Dimension]; // Serialized k dimxdim matrices in a grid

	WeightsMatrix   m_Phi;
	WeightsMatrix   m_Phi_inverse;
	WeightsMatrix   m_Phi_valid;
	WeightsMatrix   m_S;
	WeightsMatrix   m_SPrime[Dimension];

	KernelFunctionPointer m_KernelFunction;
	KernelFunctionPointer m_DerivativeKernel;
	//KernelFunctionPointer m_SecondDerivativeKernel;
	FieldPointer          m_Field;
	//FieldPointer          m_CoefficientsField;
	FieldPointer          m_GradientField;

	virtual void ComputeMatrix( WeightsMatrixType type, size_t dim = 0 );
	virtual void AfterComputeMatrix( WeightsMatrixType type );
	virtual size_t ComputeRegionOfPoint(const PointType& point, VectorType& cvector, IndexType& start, IndexType& end, OffsetTableType offsetTable );

	/** Support processing data in multiple threads. */
	itk::MultiThreader::Pointer m_Threader;
	itk::ThreadIdType           m_NumberOfThreads;

private:
	SparseMatrixTransform( const Self & );
	void operator=( const Self & );
};
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SparseMatrixTransform.hxx"
#endif
#endif /* SPARSEMATRIXTRANSFORM_H_ */
