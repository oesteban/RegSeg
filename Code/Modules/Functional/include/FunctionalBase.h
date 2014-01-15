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

#ifndef LEVELSETSBASE_H_
#define LEVELSETSBASE_H_

#include <itkObject.h>
#include <itkNumericTraits.h>
#include <itkVector.h>
#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
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

#include "CopyQuadEdgeMeshFilter.h"
#include "CopyQEMeshStructureFilter.h"
#include "WarpQEMeshFilter.h"

#include "BSplineSparseMatrixTransform.h"
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


template <typename TReferenceImageType, typename TCoordRepType = double>
class FunctionalBase: public itk::Object {
public:
	typedef FunctionalBase                    Self;
	typedef itk::Object                      Superclass;
	typedef itk::SmartPointer<Self>          Pointer;
	typedef itk::SmartPointer< const Self >  ConstPointer;

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
	typedef typename ContourType::PointDataContainerPointer  PointDataContainerPointer;
	typedef typename ContourType::PointsContainerIterator    PointsIterator;
	typedef typename ContourType::PointsContainerConstIterator    PointsConstIterator;
	typedef typename ContourType::PointIdentifier            PointIdentifier;
	typedef typename ContourType::QEType                     QEType;
	typedef typename ContourType::CellIdentifier             CellIdentifier;
	typedef typename ContourType::CellType                   CellType;
	typedef typename itk::QuadEdgeMeshPolygonCell<CellType>  PolygonType;
	typedef itk::TriangleHelper< ContourPointType >          TriangleType;

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


	typedef itk::DisplacementFieldTransform<TCoordRepType, Dimension>
															 DisplacementTransformType;
	typedef typename DisplacementTransformType::Pointer      DisplacementTransformPointer;

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

	typedef typename itk::MeshSpatialObject<ContourType>     ContourSpatialObject;
	typedef typename ContourSpatialObject::Pointer           ContourSpatialPointer;
	typedef typename ContourSpatialObject::ConstPointer      ContourSpatialConstPointer;
	typedef typename std::vector<ContourSpatialPointer>      SpatialObjectsVector;

	typedef itk::SpatialObjectToImageFilter
			       < ContourSpatialObject, ROIType >         SpatialObjectToImageFilterType;
	typedef typename SpatialObjectToImageFilterType::Pointer SpatialObjectToImageFilterPointer;

	typedef BSplineSparseMatrixTransform<TCoordRepType, Dimension, 3 > FieldInterpolatorType;
	typedef typename FieldInterpolatorType::Pointer          FieldInterpolatorPointer;

	typedef DownsampleAveragingFilter
			                 <ROIType, ProbabilityMapType >  ResampleROIFilterType;
	typedef typename ResampleROIFilterType::Pointer          ResampleROIFilterPointer;
	typedef std::vector< ResampleROIFilterPointer >          ResampleROIFilterList;


	//typedef typename itk::WarpMeshFilter
	//		< ContourType,
	//		  ContourType,
	//		  FieldType>                                     WarpContourType;

	typedef WarpQEMeshFilter< ContourType, ContourType,
			                         FieldInterpolatorType > WarpContourFilterType;
	typedef typename WarpContourFilterType::Pointer          WarpContourPointer;
	typedef std::vector< WarpContourPointer >                WarpContourFilterList;


	typedef typename std::vector< ROIPixelType >             ContourOuterRegions;
	typedef typename std::vector< ContourOuterRegions >      ContourOuterRegionsList;
	typedef typename std::vector< PointType >                PointsVector;
	typedef typename std::vector< PointsVector >             PointsList;


	typedef itk::QuadEdgeMesh< float, Dimension >            ShapeGradientType;
	typedef typename ShapeGradientType::Pointer              ShapeGradientPointer;
	typedef typename ShapeGradientType::PointDataContainer   ShapeGradientsContainer;
    typedef typename ShapeGradientsContainer::Pointer        ShapeGradientsContainerPointer;
	typedef typename std::vector<ShapeGradientPointer>       ShapeGradientList;
	typedef typename itk::CopyQEMeshStructureFilter
			                <ContourType,ShapeGradientType>  ShapeCopyType;
	typedef typename ShapeCopyType::Pointer                  ShapeCopyPointer;


	typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter< ShapeGradientType >  ContourWriterType;
	typedef typename ContourWriterType::Pointer    ContourWriterPointer;

	struct GradientSample {
		PointValueType grad;
		size_t cid;
		size_t gid;

		GradientSample() {}
		GradientSample( PointValueType g, size_t i, size_t j ): grad(g), cid(i), gid(j) {}
		GradientSample( const GradientSample &s ): grad(s.grad), cid(s.cid), gid(s.gid) {}

		GradientSample operator+(const GradientSample& g) const {
			return GradientSample( grad+g.grad, cid, gid );
		}

		GradientSample operator<(const GradientSample& g) const {
			return grad < g.grad;
		}

		GradientSample operator>(const GradientSample& g) const {
			return grad > g.grad;
		}
	};

	struct by_grad {
		bool operator()( GradientSample const &a, GradientSample const &b ) {
			return a.grad < b.grad;
		}
	};
	typedef typename std::vector< GradientSample >           SampleType;


	void CopyInformation( const FieldType* field);

	MeasureType GetValue();

	void ComputeDerivative( void );

	virtual void Initialize( void );

	ROIConstPointer GetCurrentRegion( size_t idx );
	ROIConstPointer GetCurrentRegions( void );

	const ProbabilityMapType* GetCurrentMap( size_t idx );


	size_t AddShapePrior( const ContourType* prior );
	//itkGetMacro(ShapePrior, ContourList);

	itkGetMacro(CurrentContours, ContourList);

	itkSetObjectMacro(Derivative, FieldType);
	itkGetConstObjectMacro(Derivative, FieldType);

	itkSetConstObjectMacro(ReferenceImage, ReferenceImageType);
	itkGetConstObjectMacro(ReferenceImage, ReferenceImageType);

	itkSetObjectMacro(FieldInterpolator, FieldInterpolatorType);
	itkGetObjectMacro(FieldInterpolator, FieldInterpolatorType);

protected:
	FunctionalBase();
	virtual ~FunctionalBase() {}

	void PrintSelf(std::ostream & os, itk::Indent indent) const {
		Superclass::PrintSelf(os, indent);
		os << indent << "Value: " << m_Value << std::endl;
		os << std::endl;
	}

	void InitializeSamplingGrid( void );

	inline virtual MeasureType GetEnergyOfSample( ReferencePixelType sample, size_t roi ) = 0;
	inline virtual MeasureType GetEnergyAtPoint( PointType& point, size_t roi ) = 0;
	inline virtual MeasureType GetEnergyAtPoint( PointType& point, size_t roi, ReferencePixelType& value ) = 0;


	inline bool IsInside( const PointType p, ContinuousIndex& idx ) const;
	//bool CheckExtents( const ContourType* prior ) const;

	mutable MeasureType m_Value;
	FieldPointer m_Derivative;
	FieldPointer m_ReferenceSamplingGrid;
	ContourList m_CurrentContours;
	ShapeGradientList m_Gradients;
	ConstContourList m_Priors;
	NormalFilterList m_NormalFilter;
	WarpContourPointer m_ContourUpdater;
	FieldInterpolatorPointer m_FieldInterpolator;
	DisplacementResamplerPointer m_EnergyResampler;
	ROIList m_ROIs;
	ROIList m_CurrentROIs;
	ProbabilityMapList m_CurrentMaps;
	ROIPointer m_CurrentRegions;
	DisplacementTransformPointer m_Transform;
	ReferenceImageConstPointer m_ReferenceImage;
	ContourOuterRegionsList m_OuterList;
	bool m_EnergyUpdated;
	bool m_RegionsUpdated;
	size_t m_NumberOfContours;
	size_t m_NumberOfPoints;
	size_t m_NumberOfNodes;
	double m_Scale;

	ReferencePointType m_Origin, m_End;
	DirectionType m_Direction;
	unsigned int m_SamplingFactor;

	VectorInterpolatorPointer m_LinearInterpolator;

	WarpContourFilterList m_WarpContourFilter;
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
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "FunctionalBase.hxx"
#endif

#endif /* LEVELSETSBASE_H_ */
