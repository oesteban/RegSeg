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

    using typename rstk::CachedMatrixTransform< TScalar, NDimensions >::InterpolateModeType;

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

	virtual const WeightsMatrix*  GetPhi (const bool onlyvalid = true);
	virtual const WeightsMatrix*  GetS (){
		return &this->m_S;
	}

    void SetCoefficientsImages( const CoefficientsImageArray & images );
    void SetCoefficientsImage( size_t dim, const CoefficientsImageType* c );
    void SetCoefficientsVectorImage( const FieldType* f );
    void AddCoefficientsVectorImage( const FieldType* f );
    void SetControlGridInformation( const DomainBase* image );

	void Initialize();

	void Interpolate() { this->InterpolatePoints(); this->InterpolateField(); }
	void InterpolatePoints();
	void InterpolateField();
	AltCoeffPointer GetFlatParameters();

    /** Return the multithreader used by this class. */
    itk::MultiThreader * GetMultiThreader() const { return m_Threader; }
    itkSetClampMacro( NumberOfThreads, itk::ThreadIdType, 1, ITK_MAX_THREADS);
    itkGetConstReferenceMacro(NumberOfThreads, itk::ThreadIdType);
protected:
	SparseMatrixTransform();
	~SparseMatrixTransform(){};

	enum WeightsMatrixType { PHI, PHI_FIELD, S, SPRIME, PHI_INV };

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

	CoefficientsImageArray       m_Derivatives;

	//DimensionParameters m_CoeffDerivative;  // Serialized k values in a grid
	//DimensionVector m_Jacobian[Dimension][Dimension]; // Serialized k dimxdim matrices in a grid

	WeightsMatrix   m_Phi;
	WeightsMatrix   m_Phi_inverse;
	WeightsMatrix   m_Phi_valid;
	WeightsMatrix   m_FieldPhi;
	WeightsMatrix   m_S;
	WeightsMatrix   m_SPrime[Dimension];

	KernelFunctionPointer m_KernelFunction;
	KernelFunctionPointer m_DerivativeKernel;
	//KernelFunctionPointer m_SecondDerivativeKernel;

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
