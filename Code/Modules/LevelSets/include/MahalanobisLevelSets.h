// --------------------------------------------------------------------------
// File:             MahalanobisLevelSets.h
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

#ifndef MAHALANOBISLEVELSETS_H_
#define MAHALANOBISLEVELSETS_H_

// Include headers
#include <itkObject.h>
#include <itkVector.h>
#include <itkVariableSizeMatrix.h>
#include <itkNormalQuadEdgeMeshFilter.h>
#include <itkImageToListSampleAdaptor.h>
#include <itkWeightedCovarianceSampleFilter.h>

#include "LevelSetsBase.h"


// Namespace declaration
namespace rstk {
/** \class MahalanobisLevelSets
 *  \brief This class implements a LevelSets function for registration based on
 *         Mahalanobis distance
 *
 *  MahalanobisLevelSets provides the implementation for using the Mahalanobis
 *  distance as distance function in the LevelSets environment
 *
 *  \ingroup
 */
template <typename TReferenceImageType, typename TCoordRepType = double>
class MahalanobisLevelSets: public rstk::LevelSetsBase< TReferenceImageType, TCoordRepType > {
public:
	typedef MahalanobisLevelSets                         Self;
	typedef rstk::LevelSetsBase
		< TReferenceImageType, TCoordRepType>            Superclass;
	typedef itk::SmartPointer<Self>                      Pointer;
	typedef itk::SmartPointer<const Self>                ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( MahalanobisLevelSets, rstk::LevelSetsBase );
	itkNewMacro( Self );

	typedef typename Superclass::MeasureType                 MeasureType;
	typedef typename Superclass::PointType                   PointType;
	typedef typename Superclass::PixelPointType              PixelPointType;
	typedef typename Superclass::PointValueType              PointValueType;
	typedef typename Superclass::VectorType                  VectorType;
	typedef typename Superclass::ContinuousIndex             ContinuousIndex;
	typedef typename Superclass::PixelType                   PixelType;
	typedef typename Superclass::PixelValueType              PixelValueType;
	typedef typename Superclass::DeformationFieldType        DeformationFieldType;
	typedef typename Superclass::DeformationFieldPointer     DeformationFieldPointer;
	typedef typename Superclass::ContourDeformationType      ContourDeformationType;
	typedef typename Superclass::ContourDeformationPointer   ContourDeformationPointer;
	typedef typename Superclass::ContourCopyType             ContourCopyType;
	typedef typename Superclass::ContourCopyPointer          ContourCopyPointer;
	typedef typename Superclass::NormalFilterType            NormalFilterType;
	typedef typename Superclass::NormalFilterPointer         NormalFilterPointer;

	typedef typename Superclass::ReferenceImageType          ReferenceImageType;
	typedef typename Superclass::Pointer                     ReferenceImagePointer;
	typedef typename Superclass::ConstPointer                ReferenceImageConstPointer;
	typedef typename Superclass::InterpolatorType            InterpolatorType;
	typedef typename Superclass::InterpolatorPointer         InterpolatorPointer;

	typedef typename Superclass::ROIPixelType                ROIPixelType;
	typedef typename Superclass::ROIType                     ROIType;
	typedef typename Superclass::ROIPointer                  ROIPointer;
	typedef typename Superclass::ROIConstPointer             ROIConstPointer;
	typedef typename Superclass::ROIInterpolatorType         ROIInterpolatorType;
	typedef typename Superclass::ROIList                     ROIList;

	typedef typename Superclass::ProbabilityMapType          ProbabilityMapType;
	typedef typename Superclass::ProbabilityMapPointer       ProbabilityMapPointer;
	typedef typename Superclass::ProbabilityMapConstPointer  ProbabilityMapConstPointer;
	typedef typename Superclass::ProbabilityMapList          ProbabilityMapList;

	typedef itk::Statistics::ImageToListSampleAdaptor
			< ReferenceImageType >                           ReferenceSampleType;
	typedef typename ReferenceSampleType::Pointer            ReferenceSamplePointer;

	typedef itk::Statistics::WeightedCovarianceSampleFilter
			< ReferenceSampleType >                          CovarianceFilter;
	typedef typename CovarianceFilter::Pointer               CovarianceFilterPointer;


	typedef itk::ResampleImageFilter
			                 <ROIType, ProbabilityMapType >  ResampleROIFilterType;
	typedef typename ResampleROIFilterType::Pointer          ResampleROIFilterPointer;


	itkStaticConstMacro( Components, unsigned int, itkGetStaticConstMacro(PixelType::Dimension) );
	itkStaticConstMacro( Dimension, unsigned int, TReferenceImageType::ImageDimension );

	typedef PixelType                                      MeanType;
	typedef itk::Matrix
			< PixelValueType, Components, Components >     CovarianceType;

	struct ParametersType {
		MeanType mean;
		CovarianceType cov;
		CovarianceType invcov;
		float bias;
		bool initialized;
	};

	typedef typename std::vector< ParametersType >         ParametersList;

	void SetParameters( size_t roi, ParametersType& params );

	size_t AddShapePrior( ContourDeformationType* prior, ParametersType&params );
	size_t AddShapePrior( ContourDeformationType* prior );

	void Initialize( void );

protected:
	MahalanobisLevelSets();
	~MahalanobisLevelSets() {}

	void PrintSelf( std::ostream& os, itk::Indent indent) const;

	void InitializeSamplingGrid( void );

	inline MeasureType GetEnergyAtPoint( PixelPointType& point, size_t roi );

	ParametersType UpdateParametersOfRegion( const size_t idx );
	void ComputeParameters( void );

	ParametersList m_Parameters;

private:
	MahalanobisLevelSets( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	bool CheckExtents( ContourDeformationType* prior ) const;
	bool ParametersInitialized() const;

	InterpolatorPointer m_Interp;

}; // End of class MahalanobisLevelSets
} // End of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "MahalanobisLevelSets.hxx"
#endif

#endif /* MAHALANOBISLEVELSETS_H_ */
