// --------------------------------------------------------------------------
// File:             FunctionalBase.h
// Date:             27/10/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWERegistration-Debug@Debug
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef FUNCTIONALBASE_H_
#define FUNCTIONALBASE_H_

#include <itkObject.h>
#include <itkNumericTraits.h>
#include <itkVector.h>
#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkQuadEdgeMeshTraits.h>
#include <itkQuadEdgeMesh.h>
#include <itkNormalQuadEdgeMeshFilter.h>
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

#include "rstkMacro.h"
#include "ConfigurableObject.h"
#include "CopyQuadEdgeMeshFilter.h"
#include "CopyQEMeshStructureFilter.h"
#include "WarpQEMeshFilter.h"

#include "SparseMatrixTransform.h"
#include "DownsampleAveragingFilter.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

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

	itkStaticConstMacro( Dimension, unsigned int, TReferenceImageType::ImageDimension );

	typedef double                                           MeasureType;
	typedef TCoordRepType                                    PointValueType;
	typedef itk::Vector< PointValueType, Dimension >         VectorType;
	typedef itk::Point< PointValueType, Dimension >          PointType;
	typedef itk::ContinuousIndex<TCoordRepType, Dimension >  ContinuousIndex;

	typedef TReferenceImageType                              ReferenceImageType;
	typedef typename ReferenceImageType::Pointer             ReferenceImagePointer;
	typedef typename ReferenceImageType::ConstPointer        ReferenceImageConstPointer;
	typedef typename ReferenceImageType::PixelType           ReferencePixelType;
	typedef typename ReferencePixelType::ValueType           ReferenceValueType;
	typedef typename ReferenceImageType::PointType           ReferencePointType;
	typedef typename ReferenceImageType::DirectionType       DirectionType;
	typedef typename ReferenceImageType::SizeType            ReferenceSizeType;
	typedef typename ReferenceImageType::SpacingType         ReferenceSpacingType;

	typedef itk::Array< MeasureType >                        MeasureArray;

	typedef itk::SmoothingRecursiveGaussianImageFilter< ReferenceImageType >
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 SmoothingFilterType;
	typedef typename SmoothingFilterType::Pointer			 SmoothingFilterPointer;
	typedef typename SmoothingFilterType::SigmaArrayType     SigmaArrayType;

	typedef itk::VectorLinearInterpolateImageFunction
			< ReferenceImageType >                           InterpolatorType;
	typedef typename InterpolatorType::Pointer               InterpolatorPointer;

	typedef itk::QuadEdgeMesh< VectorType, Dimension >       ContourType;
	typedef typename ContourType::Pointer                    ContourPointer;
	typedef typename ContourType::PointType                  ContourPointType;
	typedef typename ContourType::ConstPointer               ContourConstPointer;
	typedef typename std::vector<ContourPointer>             ContourList;
	typedef typename std::vector<ContourConstPointer>        ConstContourList;
	typedef typename ContourType::PointsContainerPointer     PointsContainerPointer;
	typedef typename ContourType::PointDataContainer         PointDataContainer;
	typedef typename ContourType::PointDataContainerPointer  PointDataContainerPointer;
	typedef typename ContourType::PointsContainerIterator    PointsIterator;
	typedef typename ContourType::PointsContainerConstIterator    PointsConstIterator;
	typedef typename ContourType::PointIdentifier            PointIdentifier;
	typedef typename ContourType::QEType                     QEType;
	typedef typename ContourType::CellIdentifier             CellIdentifier;
	typedef typename ContourType::CellType                   CellType;
	typedef typename itk::QuadEdgeMeshPolygonCell<CellType>  PolygonType;
	typedef itk::TriangleHelper< ContourPointType >          TriangleType;

	typedef vnl_sparse_matrix< PointValueType >              SparseMatrix;
	typedef vnl_vector< PointValueType >                     VNLVector;
	typedef itk::FixedArray< VNLVector, Dimension >          VNLVectorContainer;

	typedef itk::NormalQuadEdgeMeshFilter
			< ContourType, ContourType >                     NormalFilterType;
	typedef typename NormalFilterType::Pointer               NormalFilterPointer;
	typedef std::vector< NormalFilterPointer >               NormalFilterList;

	typedef typename itk::CopyQuadEdgeMeshFilter
			                  <ContourType,ContourType>      ContourCopyType;
	typedef typename ContourCopyType::Pointer                ContourCopyPointer;

	typedef itk::Image< VectorType, Dimension >              FieldType;
	typedef typename FieldType::Pointer                      FieldPointer;
	typedef typename FieldType::ConstPointer                 FieldConstPointer;
	typedef typename FieldType::PointType                    FieldPointType;
	typedef typename FieldType::IndexType                    FieldIndexType;

	typedef itk::VectorLinearInterpolateImageFunction
			           < FieldType >                         VectorInterpolatorType;
	typedef typename VectorInterpolatorType::Pointer         VectorInterpolatorPointer;

	typedef itk::VectorResampleImageFilter
                                  <FieldType,FieldType >     DisplacementResamplerType;
	typedef typename DisplacementResamplerType::Pointer      DisplacementResamplerPointer;

	typedef itk::DisplacementFieldJacobianDeterminantFilter
		< FieldType >                                        ModulateFilterType;
	typedef typename ModulateFilterType::Pointer             ModulateFilterPointer;

	typedef unsigned char                                    ROIPixelType;
	typedef itk::Image< ROIPixelType, Dimension >            ROIType;
	typedef typename ROIType::Pointer                        ROIPointer;
	typedef typename ROIType::ConstPointer                   ROIConstPointer;
	typedef typename itk::NearestNeighborInterpolateImageFunction
			< ROIType, TCoordRepType >                       ROIInterpolatorType;
	typedef std::vector< ROIConstPointer >                   ROIList;
	typedef itk::TriangleMeshToBinaryImageFilter
			          <ContourType, ROIType>	             BinarizeMeshFilterType;
	typedef typename BinarizeMeshFilterType::Pointer         BinarizeMeshFilterPointer;

	typedef itk::Image< float, Dimension >                   ProbabilityMapType;
	typedef typename ProbabilityMapType::Pointer             ProbabilityMapPointer;
	typedef typename ProbabilityMapType::ConstPointer        ProbabilityMapConstPointer;
	typedef std::vector< ProbabilityMapPointer >             ProbabilityMapList;

	typedef itk::NearestNeighborInterpolateImageFunction
			< ProbabilityMapType >                           MaskInterpolatorType;
	typedef typename MaskInterpolatorType::Pointer           MaskInterpolatorPointer;
	typedef itk::ResampleImageFilter
			< ProbabilityMapType, ProbabilityMapType >       ProbmapResampleType;
	typedef typename ProbmapResampleType::Pointer            ProbmapResamplePointer;

	typedef typename itk::MeshSpatialObject<ContourType>     ContourSpatialObject;
	typedef typename ContourSpatialObject::Pointer           ContourSpatialPointer;
	typedef typename ContourSpatialObject::ConstPointer      ContourSpatialConstPointer;
	typedef typename std::vector<ContourSpatialPointer>      SpatialObjectsVector;

	typedef itk::SpatialObjectToImageFilter
			       < ContourSpatialObject, ROIType >         SpatialObjectToImageFilterType;
	typedef typename SpatialObjectToImageFilterType::Pointer SpatialObjectToImageFilterPointer;

	typedef DownsampleAveragingFilter
			                 <ROIType, ProbabilityMapType >  ResampleROIFilterType;
	typedef typename ResampleROIFilterType::Pointer          ResampleROIFilterPointer;
	typedef std::vector< ResampleROIFilterPointer >          ResampleROIFilterList;

	typedef typename std::vector< ROIPixelType >             ContourOuterRegions;
	typedef typename std::vector< ContourOuterRegions >      ContourOuterRegionsList;
	typedef typename std::vector< PointType >                PointsVector;
	typedef typename std::vector< PointsVector >             PointsList;


//	typedef itk::QuadEdgeMeshTraits< float, Dimension, bool, bool >
//															 ShapeGradientTraits;
	typedef itk::QuadEdgeMesh< VectorType, Dimension> 	  	 ShapeGradientType;
	typedef typename ShapeGradientType::Pointer              ShapeGradientPointer;
	typedef typename ShapeGradientType::PointDataContainer   ShapeGradientsContainer;
    typedef typename ShapeGradientsContainer::Pointer        ShapeGradientsContainerPointer;
	typedef typename std::vector<ShapeGradientPointer>       ShapeGradientList;
	typedef typename itk::CopyQEMeshStructureFilter
			                <ContourType,ShapeGradientType>  ShapeCopyType;
	typedef typename ShapeCopyType::Pointer                  ShapeCopyPointer;

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

	itkGetMacro( CurrentContours, ContourList);
	itkGetMacro( Gradients, ShapeGradientList );
	itkGetMacro( NodesPosition, PointsVector );

	virtual void SetCurrentDisplacements( const VNLVectorContainer& vals );

	itkGetConstObjectMacro(ReferenceImage, ReferenceImageType);
	virtual void SetReferenceImage (const ReferenceImageType * _arg);

	itkGetMacro( ApplySmoothing, bool );
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
	VNLVectorContainer ComputeDerivative();
	virtual void Initialize();
	virtual void UpdateDescriptors() = 0;
	virtual std::string PrintFormattedDescriptors() = 0;

	ROIConstPointer GetCurrentRegion( size_t idx );
	itkGetConstObjectMacro( CurrentRegions, ROIType );

	itkGetConstObjectMacro( BackgroundMask, ProbabilityMapType);
	virtual void SetBackgroundMask (const ProbabilityMapType * _arg);

	itkGetConstMacro( OffMaskNodes, std::vector<size_t>);

	const ProbabilityMapType* GetCurrentMap( size_t idx );

	size_t AddShapePrior( const ContourType* prior );

	static void AddOptions( SettingsDesc& opts );

protected:
	FunctionalBase();
	virtual ~FunctionalBase() {}

	void PrintSelf(std::ostream & os, itk::Indent indent) const {
		Superclass::PrintSelf(os, indent);
		os << indent << "Value: " << m_Value << std::endl;
		os << std::endl;
	}

	void InitializeSamplingGrid( void );

	virtual MeasureType GetEnergyOfSample( ReferencePixelType sample, size_t roim, bool bias = false ) const = 0;
	MeasureType GetEnergyAtPoint( PointType& point, size_t roi ) const;
	MeasureType GetEnergyAtPoint( PointType& point, size_t roi, ReferencePixelType& value ) const;
	MeasureType EvaluateGradient( PointType& point, size_t outer_roi, size_t inner_roi ) const;


	inline bool CheckExtent( ContourPointType& p, ContinuousIndex& idx ) const;
	virtual void ParseSettings();

	size_t m_NumberOfContours;
	size_t m_NumberOfRegions;
	size_t m_NumberOfPoints;
	size_t m_NumberOfNodes;
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
	ContourList m_CurrentContours;
	ShapeGradientList m_Gradients;
	ConstContourList m_Priors;
	NormalFilterList m_NormalFilter;
	//WarpContourPointer m_ContourUpdater;
	//TransformPointer m_Transform;
	DisplacementResamplerPointer m_EnergyResampler;
	ROIList m_ROIs;
	ROIList m_CurrentROIs;
	ProbabilityMapList m_CurrentMaps;
	ProbabilityMapConstPointer m_BackgroundMask;
	ROIPointer m_CurrentRegions;
	ReferenceImageConstPointer m_ReferenceImage;
	ContourOuterRegionsList m_OuterList;
	ReferencePointType m_Origin, m_End, m_FirstPixelCenter, m_LastPixelCenter;
	ReferencePointType m_OldOrigin;
	ReferenceSizeType m_ReferenceSize;
	ReferenceSpacingType m_ReferenceSpacing;
	DirectionType m_Direction;
	DirectionType m_OldDirection;
	//VectorInterpolatorPointer m_LinearInterpolator;
	//WarpContourFilterList m_WarpContourFilter;


	InterpolatorPointer m_Interp;
	MaskInterpolatorPointer m_MaskInterp;
	PointDataContainerPointer m_CurrentDisplacements;
	PointsVector m_NodesPosition;

	std::vector<size_t> m_OffMaskNodes;
private:
	FunctionalBase(const Self &);  //purposely not implemented
	void operator=(const Self &); //purposely not implemented

	void UpdateContour();
	void ComputeCurrentRegions( void );
	void ComputeOuterRegions( void );
	void InitializeCurrentContours( void );
	void InitializeInterpolatorGrid( void );
	double ComputePointArea( const PointIdentifier &iId, ContourType *mesh );

}; // end FunctionalBase Class

itkEventMacro(WarningEvent, itk::AnyEvent);

} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "FunctionalBase.hxx"
#endif

#endif /* FUNCTIONALBASE_H_ */
