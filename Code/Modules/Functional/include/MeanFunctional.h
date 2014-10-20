// --------------------------------------------------------------------------------------
// File:          MeanFunctional.h
// Date:          Feb 10, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
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

#ifndef MEANFUNCTIONAL_H_
#define MEANFUNCTIONAL_H_

#include "FunctionalBase.h"

// Include headers
#include <itkObject.h>
#include <itkVector.h>
#include <itkVariableSizeMatrix.h>
#include <itkImageToListSampleAdaptor.h>
#include <itkWeightedCovarianceSampleFilter.h>

// Namespace declaration
namespace rstk {
/** \class MeanFunctional
 *  \brief This class implements a Functional function for registration based on
 *         Mean distance
 *
 *  MeanFunctional provides the implementation for using the Mean
 *  distance as distance function in the Functional environment
 *
 *  \ingroup
 */
template <typename TReferenceImageType, typename TCoordRepType = double>
class MeanFunctional: public rstk::FunctionalBase< TReferenceImageType, TCoordRepType > {
public:
	typedef MeanFunctional                         Self;
	typedef rstk::FunctionalBase
		< TReferenceImageType, TCoordRepType>            Superclass;
	typedef itk::SmartPointer<Self>                      Pointer;
	typedef itk::SmartPointer<const Self>                ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( MeanFunctional, rstk::FunctionalBase );
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
	MeanFunctional(){}
	~MeanFunctional(){}

	void PrintSelf( std::ostream& os, itk::Indent indent) const;

	//void InitializeSamplingGrid( void );

	inline MeasureType GetEnergyOfSample( ReferencePixelType sample, size_t roi, bool bias = false ) const;

	ParametersType UpdateParametersOfRegion( const size_t idx );
	ParametersList m_Parameters;

private:
	MeanFunctional( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	bool ParametersInitialized();

}; // End of class MeanFunctional
} // End of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "MeanFunctional.hxx"
#endif


#endif /* MEANFUNCTIONAL_H_ */
