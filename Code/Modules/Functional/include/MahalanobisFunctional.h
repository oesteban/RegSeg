// --------------------------------------------------------------------------
// File:             MahalanobisFunctional.h
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

#include "FunctionalBase.h"

// Include headers
#include <itkObject.h>
#include <itkVector.h>
#include <itkVariableSizeMatrix.h>
#include <itkImageToListSampleAdaptor.h>
#include <itkWeightedCovarianceSampleFilter.h>

// Namespace declaration
namespace rstk {
/** \class MahalanobisFunctional
 *  \brief This class implements a Functional function for registration based on
 *         Mahalanobis distance
 *
 *  MahalanobisFunctional provides the implementation for using the Mahalanobis
 *  distance as distance function in the Functional environment
 *
 *  \ingroup
 */
template <typename TReferenceImageType, typename TCoordRepType = float >
class MahalanobisFunctional: public rstk::FunctionalBase< TReferenceImageType, TCoordRepType > {
public:
	typedef MahalanobisFunctional                         Self;
	typedef rstk::FunctionalBase
		< TReferenceImageType, TCoordRepType>            Superclass;
	typedef itk::SmartPointer<Self>                      Pointer;
	typedef itk::SmartPointer<const Self>                ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( MahalanobisFunctional, rstk::FunctionalBase );
	itkNewMacro( Self );

	typedef typename Superclass::MeasureType                 MeasureType;
	typedef typename Superclass::PointType                   PointType;
	typedef typename Superclass::VectorType                  VectorType;
	typedef typename Superclass::PointValueType              PointValueType;
	typedef typename Superclass::ContinuousIndex             ContinuousIndex;

	typedef typename Superclass::ContourType                 ContourType;
	typedef typename Superclass::FieldType        FieldType;
	typedef typename Superclass::FieldPointer     FieldPointer;
	typedef typename Superclass::ContourCopyType             ContourCopyType;
	typedef typename Superclass::ContourCopyPointer          ContourCopyPointer;
	typedef typename Superclass::NormalFilterType            NormalFilterType;
	typedef typename Superclass::NormalFilterPointer         NormalFilterPointer;


	typedef typename Superclass::ReferenceImageType          ReferenceImageType;
	typedef typename Superclass::Pointer                     ReferenceImagePointer;
	typedef typename Superclass::ReferencePixelType          ReferencePixelType;
	typedef typename Superclass::ReferencePointType          ReferencePointType;
	typedef typename Superclass::ReferenceValueType          ReferenceValueType;


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


	typedef typename Superclass::ResampleROIFilterType       ResampleROIFilterType;
	typedef typename Superclass::ResampleROIFilterPointer    ResampleROIFilterPointer;


	itkStaticConstMacro( Components, unsigned int, itkGetStaticConstMacro(ReferencePixelType::Dimension) );
	itkStaticConstMacro( Dimension, unsigned int, TReferenceImageType::ImageDimension );

	typedef ReferencePixelType                               MeanType;
	typedef itk::Matrix
	      < ReferenceValueType, Components, Components >     CovarianceType;

	struct ParametersType {
		MeanType mean;
		CovarianceType cov;
		CovarianceType invcov;
		float bias;
		bool initialized;
	};

	typedef typename std::vector< ParametersType >         ParametersList;

	void SetParameters( size_t roi, ParametersType& params );

	size_t AddShapePrior( ContourType* prior, ParametersType&params );
	size_t AddShapePrior( const ContourType* prior ) { return Superclass::AddShapePrior( prior ); }

	void Initialize();
	void UpdateDescriptors();

	std::string PrintFormattedDescriptors();
protected:
	MahalanobisFunctional():m_UseRobustEstimators(true), Superclass(){}
	~MahalanobisFunctional(){}
	void PrintSelf( std::ostream& os, itk::Indent indent) const;

	inline MeasureType GetEnergyOfSample( ReferencePixelType sample, size_t roi, bool bias ) const;
	ParametersType UpdateParametersOfRegion( const size_t idx );
	ParametersType UpdateParametersOfRegionGeneral( const size_t idx );
	ParametersType UpdateParametersOfRegionRobust( const size_t idx );
	ParametersList m_Parameters;

private:
	MahalanobisFunctional( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	bool ParametersInitialized();

	bool m_UseRobustEstimators;

}; // End of class MahalanobisFunctional
} // End of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "MahalanobisFunctional.hxx"
#endif

#endif /* MAHALANOBISLEVELSETS_H_ */
