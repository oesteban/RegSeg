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

#ifndef _MAHALANOBISDISTANCEMODEL_H_
#define _MAHALANOBISDISTANCEMODEL_H_


#include <itkImageToListSampleAdaptor.h>
#include "ModelBase.h"
#include "MahalanobisDistanceMembershipFunction.h"
#include "UniformMembershipFunction.h"
#include "WeightedCovarianceSampleFilter.h"

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

	typedef itk::WeightedCovarianceSampleFilter< ReferenceSampleType >         CovarianceFilter;
	typedef typename CovarianceFilter::Pointer                                 CovarianceFilterPointer;
	typedef typename CovarianceFilter::WeightArrayType                         WeightArrayType;

	typedef std::vector< MeasurementVectorType >                               MeansContainer;
	typedef std::vector< CovarianceMatrixType >                                CovariancesContainer;

	itkGetConstMacro(RegionOffsetContainer, MeasureTypeContainer);

	std::string PrintFormattedDescriptors() override;
	virtual void ReadDescriptorsFromFile(std::string filename) override;

	inline double Evaluate(const MeasurementVectorType & x, const RegionIdentifier roi) const override {
		if( x == m_InvalidValue )
			return 0.0;
		return this->m_Memberships[roi]->Evaluate(x);
	}

protected:
	MahalanobisDistanceModel();
	virtual ~MahalanobisDistanceModel() {}
	virtual void PrintSelf(std::ostream & os, itk::Indent indent) const override;

	void GenerateData() override;
	void PostGenerateData();
	//void BeforeThreadedGenerateData();
	//void ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId);
	//void AfterThreadedGenerateData();

	virtual void InitializeMemberships() override;

	MembershipFunctionType* GetNewFunction() override {
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
