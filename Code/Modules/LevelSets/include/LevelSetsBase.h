// --------------------------------------------------------------------------
// File:             LevelSetsBase.h
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
#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
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

//#include "SparseToDenseFieldResampleFilter.h"

#include "VectorIDWBasisFunction.h"
#include "SparseMatrixTransform.h"
#include "DownsampleAveragingFilter.h"

namespace rstk {
/** \class LevelSetsBase
 *  \brief This class is a base for LevelSets definitions in the RSTK framework
 *
 *  This class is meant to provide an interface for LevelSets definitions. Given
 *  that RSTK is a registration framework, it provides the necessary elements to
 *  use the LevelSets function as a metric for non-linear registration.
 *
 *  LevelSetsBase requires the instantiation with a target image (DWI), a moving contour
 *  (surface) from T1 space and deformation fields definitions for fixed and moving
 *  objects. This implies to have a dense deformation field for DWI space and a sparse
 *  deformation field for the contour.
 *
 *  Derived classes must provide implementation for
 *  GetValue
 *  GetLevelSetsMap
 *
 *  \ingroup LevelSets
 *  \ingroup RSTK
 */


template <typename TReferenceImageType, typename TCoordRepType = double>
class LevelSetsBase: public itk::Object {
public:
	typedef LevelSetsBase                    Self;
	typedef itk::Object                      Superclass;
	typedef itk::SmartPointer<Self>          Pointer;
	typedef itk::SmartPointer< const Self >  ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro(LevelSetsBase, itk::Object);

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

	typedef itk::VectorLinearInterpolateImageFunction
			< ReferenceImageType >                           InterpolatorType;
	typedef typename InterpolatorType::Pointer               InterpolatorPointer;

	typedef itk::QuadEdgeMesh< VectorType, Dimension >       ContourType;
	typedef typename ContourType::Pointer                    ContourPointer;
	typedef typename ContourType::PointType                  ContourPointType;
	typedef typename ContourType::ConstPointer               ContourConstPointer;
	typedef typename std::vector<ContourPointer>             ContourList;
	typedef typename ContourType::PointDataContainerPointer  PointDataContainerPointer;

	typedef itk::NormalQuadEdgeMeshFilter
			< ContourType, ContourType >                     NormalFilterType;
	typedef typename NormalFilterType::Pointer               NormalFilterPointer;
	typedef std::vector< NormalFilterPointer >               NormalFilterList;

	typedef typename itk::QuadEdgeMeshToQuadEdgeMeshFilter
			<ContourType,ContourType>  ContourCopyType;
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
		< FieldType >                             ModulateFilterType;
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

	typedef typename itk::MeshSpatialObject
			                   <ContourType>      ContourSpatialObject;
	typedef typename ContourSpatialObject::Pointer           ContourSpatialPointer;
	typedef typename ContourSpatialObject::ConstPointer      ContourSpatialConstPointer;
	typedef typename std::vector<ContourSpatialPointer>      SpatialObjectsVector;

	typedef itk::SpatialObjectToImageFilter
			       < ContourSpatialObject, ROIType >         SpatialObjectToImageFilterType;
	typedef typename SpatialObjectToImageFilterType::Pointer SpatialObjectToImageFilterPointer;


	typedef RBF::VectorIDWBasisFunction< PointType, TCoordRepType, Dimension > RBFType;
	typedef SparseMatrixTransform<TCoordRepType, Dimension > FieldInterpolatorType;
	typedef typename FieldInterpolatorType::Pointer FieldInterpolatorPointer;

	typedef DownsampleAveragingFilter
			                 <ROIType, ProbabilityMapType >  ResampleROIFilterType;
	typedef typename ResampleROIFilterType::Pointer          ResampleROIFilterPointer;
	typedef std::vector< ResampleROIFilterPointer >          ResampleROIFilterList;

	//typedef SparseToDenseFieldResampleFilter
	//	<ContourType, FieldType>       SparseToDenseFieldResampleType;
	//typedef typename SparseToDenseFieldResampleType::Pointer SparseToDenseFieldResamplePointer;

	typedef typename itk::WarpMeshFilter
			< ContourType,
			  ContourType,
			  FieldType>                                     WarpContourType;
	typedef typename WarpContourType::Pointer                WarpContourPointer;

	typedef typename std::vector< ROIPixelType >             ContourOuterRegions;
	typedef typename std::vector< ContourOuterRegions >      ContourOuterRegionsList;
	typedef typename std::vector< PointType >                ControlPointsVector;
	typedef typename std::vector< ControlPointsVector >      ControlPointsList;

	void CopyInformation( const FieldType* field);

	MeasureType GetValue();

	void ComputeGradient( void );

	MeasureType UpdateContour( const FieldType* newField );

	virtual void Initialize( void );

	ROIConstPointer GetCurrentRegion( size_t idx );
	ROIConstPointer GetCurrentRegions( void );

	const ProbabilityMapType* GetCurrentMap( size_t idx );

	size_t AddShapePrior( ContourType* prior );
	//itkGetMacro(ShapePrior, ContourList);

	itkGetMacro(CurrentContourPosition, ContourList);

	itkSetObjectMacro(GradientMap, FieldType);
	itkGetConstObjectMacro(GradientMap, FieldType);

	itkSetConstObjectMacro(ReferenceImage, ReferenceImageType);
	itkGetConstObjectMacro(ReferenceImage, ReferenceImageType);

protected:
	LevelSetsBase();
	virtual ~LevelSetsBase() {}

	void PrintSelf(std::ostream & os, itk::Indent indent) const {
		Superclass::PrintSelf(os, indent);
		os << indent << "Value: " << m_Value << std::endl;
		os << std::endl;
	}

	void InitializeSamplingGrid( void );
	inline virtual MeasureType GetEnergyAtPoint( PointType& point, size_t roi ) = 0;


	inline bool IsInside( const PointType p, ContinuousIndex& idx ) const;
	//bool CheckExtents( const ContourType* prior ) const;

	mutable MeasureType m_Value;
	FieldPointer m_GradientMap;
	FieldPointer m_ReferenceSamplingGrid;
	ContourList m_CurrentContourPosition;
	NormalFilterList m_NormalFilter;
	ContourCopyPointer m_ContourCopier;
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
	ControlPointsList m_ShapePrior;
	bool m_Modified;
	bool m_RegionsModified;
	size_t m_NumberOfContours;

private:
	LevelSetsBase(const Self &);  //purposely not implemented
	void operator=(const Self &); //purposely not implemented

	void ComputeCurrentRegions( void );


}; // end LevelSetsBase Class
} // end namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "LevelSetsBase.hxx"
#endif

#endif /* LEVELSETSBASE_H_ */
