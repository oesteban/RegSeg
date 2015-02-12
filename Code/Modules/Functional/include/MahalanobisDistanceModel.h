// --------------------------------------------------------------------------------------
// File:          MahalanobisDistanceModel.h
// Date:          Dec 23, 2014
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
// This file is part of ACWEReg
//
// ACWEReg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWEReg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWEReg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef _MAHALANOBISDISTANCEMODEL_H_
#define _MAHALANOBISDISTANCEMODEL_H_


#include <itkImageToListSampleAdaptor.h>
#include "ModelBase.h"
#include "MahalanobisDistanceMembershipFunction.h"
#include "UniformMembershipFunction.h"
#include "WeightedCovarianceSampleFilter.h"

#include <itkWeightedCovarianceSampleFilter.h>

namespace rstk {

template< typename TInputVectorImage, typename TPriorsPrecisionType = float>
class MahalanobisDistanceModel: public ModelBase< TInputVectorImage, TPriorsPrecisionType > {

public:
	/** Standard class typedefs */
	typedef MahalanobisDistanceModel                                  Self;
	typedef ModelBase< TInputVectorImage, TPriorsPrecisionType >      Superclass;
	typedef itk::SmartPointer< Self >                                 Pointer;
	typedef itk::SmartPointer< const Self >                           ConstPointer;

	/** Standard macros */
	itkTypeMacro(MahalanobisDistanceModel, ModelBase);
	itkNewMacro(MahalanobisDistanceModel);

	itkStaticConstMacro(Dimension, unsigned int, Superclass::Dimension);

	typedef typename Superclass::InputImageType                                InputImageType;
	typedef typename Superclass::Pointer                                       InputImagePointer;
	typedef typename Superclass::SizeType                                      SizeType;
	typedef typename Superclass::SpacingType                                   SpacingType;
	typedef typename Superclass::IndexType                                     IndexType;
	typedef typename Superclass::PixelType                                     PixelType;
	typedef typename Superclass::PointType                                     PointType;
	typedef typename Superclass::RegionType                                    RegionType;
	typedef typename Superclass::PixelValueType                                PixelValueType;

	/** MeasurementVector typedef support */
	typedef typename Superclass::MeasureType                                   MeasureType;
	typedef typename Superclass::MeasurementVectorType                         MeasurementVectorType;
	typedef typename Superclass::MeasurementVectorSizeType                     MeasurementVectorSizeType;
	typedef typename Superclass::RegionIdentifier                              RegionIdentifier;
	typedef typename Superclass::MeasureTypeContainer                          MeasureTypeContainer;

	typedef typename Superclass::PriorsPrecisionType                           PriorsPrecisionType;
	typedef typename Superclass::PriorsImageType                               PriorsImageType;
	typedef typename Superclass::PriorsPixelType                               PriorsPixelType;
	typedef typename Superclass::PriorsImagePointer                            PriorsImagePointer;
	typedef typename Superclass::PriorsImageIteratorType                       PriorsImageIteratorType;
	typedef typename Superclass::PriorsAdaptor                                 PriorsAdaptor;
	typedef typename Superclass::PriorsAdaptorPointer                          PriorsAdaptorPointer;

	typedef typename Superclass::MaskType                                      MaskType;
	typedef typename Superclass::MaskPointer                                   MaskPointer;

	typedef typename Superclass::MembershipFunctionType                        MembershipFunctionType;
	typedef typename Superclass::MembershipFunctionPointer                     MembershipFunctionPointer;
	typedef typename Superclass::MembershipFunctionsArray                      MembershipFunctionsArray;

	typedef MahalanobisDistanceMembershipFunction<MeasurementVectorType>       InternalFunctionType;
	typedef typename InternalFunctionType::Pointer                             InternalFunctionPointer;
	typedef typename InternalFunctionType::CovarianceMatrixType                CovarianceMatrixType;

	typedef UniformMembershipFunction<MeasurementVectorType>                   UniformFunctionType;
	typedef typename UniformFunctionType::Pointer                              UniformFunctionPointer;

	typedef itk::Statistics::ImageToListSampleAdaptor< InputImageType >        ReferenceSampleType;
	typedef typename ReferenceSampleType::Pointer                              ReferenceSamplePointer;

	typedef rstk::WeightedCovarianceSampleFilter< ReferenceSampleType >        CovarianceFilter;
	typedef typename CovarianceFilter::Pointer                                 CovarianceFilterPointer;
	typedef typename CovarianceFilter::WeightArrayType                         WeightArrayType;

	typedef std::vector< MeasurementVectorType >                               MeansContainer;
	typedef std::vector< CovarianceMatrixType >                                CovariancesContainer;

	itkGetConstMacro(RegionOffsetContainer, MeasureTypeContainer);

	std::string PrintFormattedDescriptors();
	virtual void ReadDescriptorsFromFile(std::string filename);

	inline double Evaluate(const MeasurementVectorType & x, const RegionIdentifier roi) const {
		if( x == m_InvalidValue )
			return 0.0;
		return this->m_Memberships[roi]->Evaluate(x);
	}

protected:
	MahalanobisDistanceModel();
	virtual ~MahalanobisDistanceModel() {}
	virtual void PrintSelf(std::ostream & os, itk::Indent indent) const;

	void GenerateData();
	void PostGenerateData();
	//void BeforeThreadedGenerateData();
	//void ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId);
	//void AfterThreadedGenerateData();

	virtual void InitializeMemberships();

	MembershipFunctionType* GetNewFunction() {
		return dynamic_cast<MembershipFunctionType*>(InternalFunctionType::New().GetPointer());
	}

private:
	MeasureType ComputeCovarianceDeterminant(CovarianceMatrixType& cov) const;

	void Estimate();
	void EstimateRobust();

	MahalanobisDistanceModel(const Self &);   //purposely not implemented
	void operator=(const Self &);             //purposely not implemented

	MeansContainer        m_Means;
	MeansContainer        m_RangeLower;
	MeansContainer        m_RangeUpper;
	CovariancesContainer  m_Covariances;
	MeasureTypeContainer  m_RegionOffsetContainer;
	MeasurementVectorType m_InvalidValue;
};
}


#ifndef ITK_MANUAL_INSTANTIATION
#include "MahalanobisDistanceModel.hxx"
#endif

#endif /* _MAHALANOBISDISTANCEMODEL_H_ */
