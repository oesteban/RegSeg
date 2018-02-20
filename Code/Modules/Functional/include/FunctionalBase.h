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

#ifndef FUNCTIONALBASE_H_
#define FUNCTIONALBASE_H_

#include <mutex>

#include <itkObject.h>
#include <itkNumericTraits.h>
#include <itkVector.h>
#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include "VectorLinearInterpolateImageFunction.h"


#include <itkImageFileReader.h>
#include <itkMeshFileReader.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkQuadEdgeMeshTraits.h>
#include <itkQuadEdgeMesh.h>
#include <itkWarpMeshFilter.h>
#include <itkMeshSpatialObject.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <itkDisplacementFieldTransform.h>
#include <itkQuadEdgeMeshPolygonCell.h>
#include <itkTriangleHelper.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <vnl/vnl_sparse_matrix.h>

#include <itkMahalanobisDistanceMembershipFunction.h>

#include <itkSimplexMesh.h>
#include <itkSimplexMeshToTriangleMeshFilter.h>
#include <itkTriangleMeshToSimplexMeshFilter.h>


#include "rstkMacro.h"
#include "NormalQuadEdgeMeshFilter.h"
#include "ConfigurableObject.h"
#include "CopyQuadEdgeMeshFilter.h"
#include "CopyCastMeshFilter.h"
#include "CopyQEMeshStructureFilter.h"
#include "WarpQEMeshFilter.h"
#include "SparseMatrixTransform.h"
#include "DownsampleAveragingFilter.h"
#include "MultilabelBinarizeMeshFilter.h"

#include "EnergyCalculatorFilter.h"
#include "MahalanobisDistanceModel.h"

namespace rstk {
/** \class FunctionalBase
 *  \brief This class is a base for Functional definitions in the RSTK framework
 *
 *  This class is meant to provide an interface for Functional definitions. Given
 *  that RSTK is a registration framework, it provides the necessary elements to
 *  use the Functional function as a metric for non-linear registration.
 *
 *  FunctionalBase requires the instantiation with a target image (DWI), a moving contour
 *  (surface) from T1 space and deformation fields definitions for fixed and moving
 *  objects. This implies to have a dense deformation field for DWI space and a sparse
 *  deformation field for the contour.
 *
 *  Derived classes must provide implementation for
 *  GetValue
 *  GetFunctionalMap
 *
 *  \ingroup Functional
 *  \ingroup RSTK
 */


template <typename TReferenceImageType, typename TCoordRepType = float>
class FunctionalBase:
		public itk::Object,
		public ConfigurableObject  {
public:
	typedef FunctionalBase                                   Self;
	typedef itk::Object                                      Superclass;
	typedef itk::SmartPointer<Self>                          Pointer;
	typedef itk::SmartPointer< const Self >                  ConstPointer;
	typedef ConfigurableObject                               SettingsClass;
	typedef typename SettingsClass::SettingsMap              SettingsMap;
	typedef typename SettingsClass::SettingsDesc             SettingsDesc;
	typedef std::vector< SettingsMap >                       SettingsList;

	/** Run-time type information (and related methods). */
	itkTypeMacro(FunctionalBase, itk::Object);
	itkNewMacro(FunctionalBase);

	itkStaticConstMacro( Dimension, unsigned int, TReferenceImageType::ImageDimension );

	typedef double                                                    MeasureType;
	typedef TCoordRepType                                             PointValueType;
	typedef itk::Vector< PointValueType, Dimension >                  VectorType;
	typedef itk::Point< PointValueType, Dimension >                   PointType;
	typedef itk::ContinuousIndex<TCoordRepType, Dimension >           ContinuousIndex;

	typedef itk::FixedArray<TCoordRepType, Dimension >                ScalesType;

	typedef TReferenceImageType                                       ReferenceImageType;
	typedef typename ReferenceImageType::Pointer                      ReferenceImagePointer;
	typedef typename ReferenceImageType::ConstPointer                 ReferenceImageConstPointer;
	typedef typename ReferenceImageType::PixelType                    ReferencePixelType;
	typedef typename ReferenceImageType::InternalPixelType            ChannelPixelType;
	typedef typename ReferencePixelType::ValueType                    ReferenceValueType;
	typedef typename ReferenceImageType::PointType                    ReferencePointType;
	typedef typename ReferenceImageType::IndexType                    ReferenceIndexType;
	typedef typename ReferenceImageType::DirectionType                DirectionType;
	typedef typename ReferenceImageType::SizeType                     ReferenceSizeType;
	typedef typename ReferenceImageType::SpacingType                  ReferenceSpacingType;
	typedef itk::Array< MeasureType >                                 MeasureArray;

	typedef itk::SmoothingRecursiveGaussianImageFilter
			< ReferenceImageType >                                    SmoothingFilterType;
	typedef typename SmoothingFilterType::Pointer			          SmoothingFilterPointer;
	typedef typename SmoothingFilterType::SigmaArrayType              SigmaArrayType;

	typedef rstk::VectorLinearInterpolateImageFunction
			< ReferenceImageType >                                    InterpolatorType;
	typedef typename InterpolatorType::Pointer                        InterpolatorPointer;

	typedef itk::QuadEdgeMesh< VectorType, Dimension >                VectorContourType;
	typedef typename VectorContourType::Pointer                       VectorContourPointer;
	typedef typename VectorContourType::PointType                     VectorContourPointType;
	typedef typename VectorContourType::ConstPointer                  VectorContourConstPointer;
	typedef typename std::vector<VectorContourPointer>                VectorContourList;
	typedef typename std::vector<VectorContourConstPointer>           VectorConstContourList;
	typedef typename VectorContourType::PointsContainerPointer        PointsContainerPointer;
	typedef typename VectorContourType::PointDataContainer            PointDataContainer;
	typedef typename VectorContourType::PointDataContainerPointer     PointDataContainerPointer;
	typedef typename VectorContourType::PointsContainerIterator       PointsIterator;
	typedef typename VectorContourType::PointsContainerConstIterator  PointsConstIterator;
	typedef typename VectorContourType::PointIdentifier               PointIdentifier;
	typedef typename VectorContourType::QEType                        QEType;
	typedef typename VectorContourType::CellIdentifier                CellIdentifier;
	typedef typename VectorContourType::CellType                      CellType;
	typedef typename itk::QuadEdgeMeshPolygonCell<CellType>           PolygonType;
	typedef itk::TriangleHelper< VectorContourPointType >             TriangleType;

	typedef std::vector< PointIdentifier >                            PointIdContainer;
	typedef typename PointIdContainer::iterator                       PointIdIterator;

	typedef itk::QuadEdgeMesh< PointValueType, Dimension >            ScalarContourType;
	typedef typename ScalarContourType::Pointer                       ScalarContourPointer;
	typedef typename ScalarContourType::PointType                     ScalarContourPointType;
	typedef typename ScalarContourType::ConstPointer                  ScalarContourConstPointer;
	typedef typename std::vector<ScalarContourPointer>                ScalarContourList;
	typedef typename std::vector<ScalarContourConstPointer>           ScalarConstContourList;

	typedef typename itk::Image<ChannelPixelType, Dimension>          ChannelType;
	typedef typename itk::ImageFileReader<ChannelType>                ChannelReader;
	typedef typename itk::MeshFileReader<ScalarContourType>           PriorReader;

	typedef vnl_sparse_matrix< PointValueType >                       SparseMatrix;
	typedef vnl_vector< PointValueType >                              VNLVector;
	typedef itk::FixedArray< VNLVector, Dimension >                   VNLVectorContainer;

	typedef rstk::NormalQuadEdgeMeshFilter
			< VectorContourType, VectorContourType >                  NormalFilterType;
	typedef typename NormalFilterType::Pointer                        NormalFilterPointer;
	typedef typename NormalFilterType::AreaContainerType              NormalFilterAreasContainer;
	typedef std::vector< NormalFilterPointer >                        NormalFilterList;

	// Contour copiers
	typedef typename rstk::CopyQuadEdgeMeshFilter
				  <ScalarContourType, ScalarContourType>              ScalarContourCopyType;
	typedef typename ScalarContourCopyType::Pointer                   ScalarContourCopyPointer;

	typedef typename rstk::CopyQuadEdgeMeshFilter
				  <VectorContourType, VectorContourType>              VectorContourCopyType;
	typedef typename VectorContourCopyType::Pointer                   VectorContourCopyPointer;
	typedef typename rstk::CopyCastMeshFilter
			      <ScalarContourType, VectorContourType>              Scalar2VectorCopyType;
	typedef typename Scalar2VectorCopyType::Pointer                   Scalar2VectorCopyPointer;
	typedef typename rstk::CopyCastMeshFilter
			      <VectorContourType,ScalarContourType>               Vector2ScalarCopyType;
	typedef typename Vector2ScalarCopyType::Pointer                   Vector2ScalarCopyPointer;

	typedef itk::Image< VectorType, Dimension >                       FieldType;
	typedef typename FieldType::Pointer                               FieldPointer;
	typedef typename FieldType::ConstPointer                          FieldConstPointer;
	typedef typename FieldType::PointType                             FieldPointType;
	typedef typename FieldType::IndexType                             FieldIndexType;

	typedef itk::VectorLinearInterpolateImageFunction
			           < FieldType >                                  VectorInterpolatorType;
	typedef typename VectorInterpolatorType::Pointer                  VectorInterpolatorPointer;

	typedef itk::VectorResampleImageFilter
                                  <FieldType,FieldType >              DisplacementResamplerType;
	typedef typename DisplacementResamplerType::Pointer               DisplacementResamplerPointer;

	typedef itk::DisplacementFieldJacobianDeterminantFilter
		< FieldType >                                                 ModulateFilterType;
	typedef typename ModulateFilterType::Pointer                      ModulateFilterPointer;

	typedef unsigned char                                             ROIPixelType;
	typedef itk::Image< ROIPixelType, Dimension >                     ROIType;
	typedef typename ROIType::Pointer                                 ROIPointer;
	typedef typename ROIType::ConstPointer                            ROIConstPointer;
	typedef typename itk::NearestNeighborInterpolateImageFunction
			< ROIType, TCoordRepType >                                ROIInterpolatorType;
	typedef itk::ResampleImageFilter
			< ROIType, ROIType, TCoordRepType >                       ROIResampleType;
	typedef std::vector< ROIConstPointer >                            ROIList;

	typedef MahalanobisDistanceModel< ReferenceImageType >            EnergyModelType;
	typedef typename EnergyModelType::Pointer                         EnergyModelPointer;

	typedef EnergyCalculatorFilter<ReferenceImageType>                EnergyFilter;
	typedef typename EnergyFilter::Pointer                            EnergyFilterPointer;
	typedef typename EnergyFilter::PriorsImageType                    PriorsImageType;
	typedef typename PriorsImageType::Pointer                         PriorsImagePointer;
	typedef typename PriorsImageType::PixelType                       PriorsPixelType;
	typedef typename PriorsImageType::InternalPixelType               PriorsValueType;

	typedef MultilabelBinarizeMeshFilter< VectorContourType >         BinarizeMeshFilterType;
	typedef typename BinarizeMeshFilterType::Pointer                  BinarizeMeshFilterPointer;
	typedef typename BinarizeMeshFilterType::OutputImageType          BinarizationImageType;

	typedef DownsampleAveragingFilter
			< BinarizationImageType, PriorsImageType >                DownsampleFilter;
	typedef typename DownsampleFilter::Pointer                        DownsamplePointer;

	typedef itk::Image< float, Dimension >                            ProbabilityMapType;
	typedef typename ProbabilityMapType::Pointer                      ProbabilityMapPointer;
	typedef typename ProbabilityMapType::ConstPointer                 ProbabilityMapConstPointer;
	typedef std::vector< ProbabilityMapPointer >                      ProbabilityMapList;

	typedef itk::NearestNeighborInterpolateImageFunction
			< ProbabilityMapType >                                    MaskInterpolatorType;
	typedef typename MaskInterpolatorType::Pointer                    MaskInterpolatorPointer;
	typedef itk::ResampleImageFilter
			< ProbabilityMapType, ProbabilityMapType >                ProbmapResampleType;
	typedef typename ProbmapResampleType::Pointer                     ProbmapResamplePointer;

	typedef typename itk::MeshSpatialObject<VectorContourType>        ContourSpatialObject;
	typedef typename ContourSpatialObject::Pointer                    ContourSpatialPointer;
	typedef typename ContourSpatialObject::ConstPointer               ContourSpatialConstPointer;
	typedef typename std::vector<ContourSpatialPointer>               SpatialObjectsVector;

	typedef itk::SpatialObjectToImageFilter
			       < ContourSpatialObject, ROIType >                  SpatialObjectToImageFilterType;
	typedef typename SpatialObjectToImageFilterType::Pointer          SpatialObjectToImageFilterPointer;

	typedef DownsampleAveragingFilter
			                 <ROIType, ProbabilityMapType >           ResampleROIFilterType;
	typedef typename ResampleROIFilterType::Pointer                   ResampleROIFilterPointer;
	typedef std::vector< ResampleROIFilterPointer >                   ResampleROIFilterList;

	typedef typename std::vector< ROIPixelType >                      ContourOuterRegions;
	typedef typename std::vector< ContourOuterRegions >               ContourOuterRegionsList;
	typedef typename std::vector< PointType >                         PointsVector;
	typedef typename std::vector< PointsVector >                      PointsList;
	typedef typename std::vector< PointValueType >                    PointValuesVector;

	// typedef itk::QuadEdgeMesh< VectorType, Dimension> 	  	          ShapeGradientType;
	// typedef typename ShapeGradientType::Pointer                       ShapeGradientPointer;
	// typedef typename ShapeGradientType::PointDataContainer            ShapeGradientsContainer;
    // typedef typename ShapeGradientsContainer::Pointer                 ShapeGradientsContainerPointer;
	// typedef typename std::vector<ShapeGradientPointer>                ShapeGradientList;
	// typedef typename itk::CopyQEMeshStructureFilter
	// 		                <VectorContourType,ShapeGradientType>     ShapeCopyType;
	// typedef typename ShapeCopyType::Pointer                           ShapeCopyPointer;
	typedef itk::FixedArray< PointValueType, 7u >                     GradientStatsArray;

	struct GradientSample {
		PointValueType grad;
		PointValueType w;
		VectorType normal;
		size_t cid;         // point id in contour sid
		size_t gid;         // global point id
		size_t sid;         // shape id

		//GradientSample(): grad(0), cid(0), gid(0) {}
		GradientSample( PointValueType g, VectorType n, size_t i ): grad(g), w(1.0), normal(n), cid(0),gid(i),sid(0) {}
		GradientSample( PointValueType g, PointValueType weight, VectorType n, size_t i, size_t j, size_t k ): grad(g), w(weight), normal(n), cid(i), gid(j), sid(k) {}
		GradientSample( const GradientSample &s ): grad(s.grad), w(s.w), normal(s.normal), cid(s.cid), gid(s.gid), sid(s.sid) {}

		GradientSample operator+(const GradientSample& g) const {
			PointValueType val = (grad*w + g.grad * g.w )/ (w+g.w);
			return GradientSample( grad+g.grad, cid, gid );
		}

		bool operator<(const GradientSample& g) const {
			return grad < g.grad;
		}

		bool operator>(const GradientSample& g) const {
			return grad > g.grad;
		}
	};

//	struct by_grad {
//		bool operator()( GradientSample const &a, GradientSample const &b ) {
//			return a.grad < b.grad;
//		}
//	};
	struct by_gid {
		bool operator()( GradientSample const &a, GradientSample const &b ) {
			return a.gid < b.gid;
		}
	};
	struct by_cid {
		bool operator()( GradientSample const &a, GradientSample const &b ) {
			return a.cid < b.cid;
		}
	};


	typedef typename std::vector< GradientSample >           SampleType;

	itkSetClampMacro( DecileThreshold, float, 0.0, 0.5 );
	itkGetMacro( DecileThreshold, float );

	itkGetMacro( CurrentContours, VectorContourList);
	itkGetMacro( Gradients, VectorContourList );
	itkGetMacro( Vertices, PointsVector );
	itkGetMacro( ValidVertices, PointIdContainer );

	virtual void SetCurrentDisplacements( const VNLVectorContainer& vals );

	itkGetConstObjectMacro(ReferenceImage, ReferenceImageType);
	itkSetConstObjectMacro(ReferenceImage, ReferenceImageType);

	void LoadReferenceImage( const std::vector<std::string> fixedImageNames );

	itkGetConstObjectMacro( CurrentMaps, PriorsImageType);

	itkGetMacro( ApplySmoothing, bool );
	itkGetMacro( UseBackground, bool );
	itkSetMacro( UseBackground, bool );

	itkGetMacro( Sigma, SigmaArrayType );
	itkSetMacro( Sigma, SigmaArrayType );

	itkGetMacro( MaxEnergy, MeasureType );

	void SetSigma( float s ) {
		itkDebugMacro( "set Sigma to " << s );
		bool modified = false;
		for( size_t i = 0; i<Dimension; i++){
			if ( this->m_Sigma[i] != s ) {
				modified = true;
				break;
			}
		}
		if ( modified ) {
			this->m_Sigma.Fill(s);
			this->m_ApplySmoothing = true;
			this->Modified();
		}
	}

	MeasureType GetValue();
	itkGetConstMacro(RegionValue, MeasureArray);
	itkGetConstMacro(GradientStatistics, GradientStatsArray);
	void ComputeDerivative(PointValueType* gradVector, ScalesType scales);

	virtual void Initialize();
	virtual void UpdateDescriptors() {
		this->m_Model->SetPriorsMap(this->m_CurrentMaps);
		this->m_Model->Update();
		this->m_MaxEnergy = this->m_Model->GetMaxEnergy();
	}

	virtual std::string PrintFormattedDescriptors() {
		return this->m_Model->PrintFormattedDescriptors();
	}

	itkGetConstObjectMacro( CurrentRegions, ROIType );

	itkGetConstObjectMacro( BackgroundMask, ProbabilityMapType);
	virtual void SetBackgroundMask (const ProbabilityMapType * _arg);

	itkGetConstMacro( OffMaskVertices, std::vector<size_t>);

	void LoadShapePriors( std::vector< std::string > movingSurfaceNames );
	size_t AddShapePrior( const ScalarContourType* prior );

	size_t AddShapeTarget( const ScalarContourType* surf ) {
		this->m_Target.push_back( surf );
	}

	static void AddOptions( SettingsDesc& opts );

	const MeasureArray GetFinalEnergy() const;
	std::string GetInfoString() {
		this->m_InfoBuffer << " } }";
		std::string result = this->m_InfoBuffer.str();
		this->m_InfoBuffer.str("");
		this->m_InfoBuffer << "{ \"info\": { ";
		return result;
	}

    /** Return the multithreader used by this class. */
    itk::MultiThreader * GetMultiThreader() const { return m_Threader; }
    itkSetClampMacro( NumberOfThreads, itk::ThreadIdType, 1, ITK_MAX_THREADS);
    itkGetConstReferenceMacro(NumberOfThreads, itk::ThreadIdType);

protected:
	FunctionalBase();
	virtual ~FunctionalBase() {}

	void PrintSelf(std::ostream & os, itk::Indent indent) const override {
		Superclass::PrintSelf(os, indent);
		os << indent << "Value: " << m_Value << std::endl;
		os << std::endl;
	}

	void InitializeSamplingGrid( void );

	//virtual MeasureType GetEnergyOfSample( ReferencePixelType sample, size_t roim, bool bias = false ) const = 0;
	MeasureType GetEnergyAtPoint( const PointType& point, size_t roi ) const;
	MeasureType GetEnergyAtPoint( const PointType& point, size_t roi, ReferencePixelType& value ) const;
	inline MeasureType EvaluateGradient( const PointType& point, size_t outer_roi, size_t inner_roi ) const;

	inline bool CheckExtent( VectorContourPointType& p, ContinuousIndex& idx ) const;
	virtual void ParseSettings() override;
	//virtual MeasureType GetEnergyOffset(size_t roi) const = 0;

	// Methods for multithreading
	struct ParallelGradientStruct {
		Self* selfptr;
		size_t total;
		PointValuesVector* gradients;
		std::mutex mutex;
		std::vector<NormalFilterAreasContainer> areas;
		std::vector<PointsContainerPointer> points;
		std::vector< double > totalAreas;
	};

	static ITK_THREAD_RETURN_TYPE ThreadedDerivativeCallback(void *arg);
	PointValuesVector ThreadedDerivativeCompute(size_t start, size_t stop,
			std::vector<PointsContainerPointer> points,
			std::vector<NormalFilterAreasContainer> areas,
			std::vector< double > totalAreas);



	size_t m_NumberOfContours;
	size_t m_NumberOfRegions;
	size_t m_NumberOfVertices;
	size_t m_SamplingFactor;
	SigmaArrayType m_Sigma;
	float m_DecileThreshold;
	bool m_DisplacementsUpdated;
	bool m_EnergyUpdated;
	bool m_RegionsUpdated;
	bool m_ApplySmoothing;
	bool m_UseBackground;

	mutable MeasureType m_Value;
	mutable MeasureArray m_RegionValue;
	mutable MeasureType m_MaxEnergy;
	FieldPointer m_ReferenceSamplingGrid;
	VectorContourList m_CurrentContours;
	VectorContourList m_Gradients;
	ScalarConstContourList m_Priors;
	ScalarConstContourList m_Target;
	EnergyModelPointer m_Model;
	EnergyFilterPointer m_EnergyCalculator;
	// ROIList m_ROIs;
	ROIList m_CurrentROIs;
	PriorsImagePointer m_CurrentMaps;
	ProbabilityMapConstPointer m_BackgroundMask;
	ROIPointer m_CurrentRegions;
	ReferenceImageConstPointer m_ReferenceImage;
	ReferencePointType m_Origin, m_End, m_FirstPixelCenter, m_LastPixelCenter;
	ReferencePointType m_ITKOrientOffset;
	ReferencePointType m_PhyExtentMin, m_PhyExtentMax;
	ReferenceSizeType m_ReferenceSize;
	ReferenceSpacingType m_ReferenceSpacing;
	DirectionType m_Direction;
	//VectorInterpolatorPointer m_LinearInterpolator;
	//WarpContourFilterList m_WarpContourFilter;


	InterpolatorPointer m_Interp;
	MaskInterpolatorPointer m_MaskInterp;
	PointDataContainerPointer m_CurrentDisplacements;
	PointsVector m_Vertices;
	PointIdContainer m_ValidVertices;
	PointIdContainer m_OuterRegion;
	PointIdContainer m_InnerRegion;
	PointIdContainer m_Offsets;

	std::vector<size_t> m_OffMaskVertices;

	GradientStatsArray m_GradientStatistics;

	mutable std::stringstream m_InfoBuffer;

private:
	FunctionalBase(const Self &);  //purposely not implemented
	void operator=(const Self &); //purposely not implemented

	void UpdateContour();
	void ComputeCurrentRegions();
	void InitializeContours();
	void InitializeInterpolatorGrid();
	double ComputePointArea( const PointIdentifier &iId, VectorContourType *mesh );

	/** Support processing data in multiple threads. */
	itk::MultiThreader::Pointer m_Threader;
	itk::ThreadIdType           m_NumberOfThreads;
}; // end FunctionalBase Class

itkEventMacro(WarningEvent, itk::AnyEvent);

} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "FunctionalBase.hxx"
#endif

#endif /* FUNCTIONALBASE_H_ */
