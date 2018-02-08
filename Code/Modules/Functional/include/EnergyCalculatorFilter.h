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

#ifndef _ENERGYCALCULATORFILTER_H_
#define _ENERGYCALCULATORFILTER_H_

#include <itkImageTransformer.h>
#include <itkImageRegionConstIterator.h>
#include "ModelBase.h"

namespace rstk {


template < typename TInputVectorImage, typename TMeasureType = double, typename TPriorsPrecisionType = float>
class EnergyCalculatorFilter: public itk::ImageTransformer< TInputVectorImage > {

public:
	typedef EnergyCalculatorFilter                                            Self;
	typedef itk::ImageTransformer< TInputVectorImage >                        Superclass;
	typedef itk::SmartPointer< Self >                                         Pointer;
	typedef itk::SmartPointer< const Self >                                   ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(EnergyCalculatorFilter, itk::ImageTransformer);
	itkStaticConstMacro(Dimension, unsigned int, TInputVectorImage::ImageDimension);

	typedef TInputVectorImage                                                 InputImageType;
	typedef typename InputImageType::Pointer                                  InputImagePointer;
	typedef typename InputImageType::SizeType                                 SizeType;
	typedef typename InputImageType::SpacingType                              SpacingType;
	typedef typename InputImageType::IndexType                                IndexType;
	typedef typename InputImageType::PixelType                                PixelType;
	typedef typename InputImageType::PointType                                PointType;
	typedef typename InputImageType::RegionType                               RegionType;

	typedef TPriorsPrecisionType                                              PriorsPrecisionType;
	typedef itk::VectorImage< PriorsPrecisionType, Dimension >                PriorsImageType;
	typedef typename PriorsImageType::PixelType                               PriorsPixelType;
	typedef typename PriorsImageType::Pointer                                 PriorsImagePointer;
	typedef itk::ImageRegionConstIterator< PriorsImageType >                  PriorsImageIteratorType;
	typedef itk::Array<PriorsPrecisionType>                                   TotalVolumeContainer;

	typedef itk::Image< PriorsPrecisionType, Dimension >                      MaskType;
	typedef typename MaskType::Pointer                                        MaskPointer;

	typedef TMeasureType                                                      MeasureType;
	typedef itk::Array<MeasureType>                                           MeasureArrayType;
	typedef itk::SimpleDataObjectDecorator< MeasureArrayType >                MeasureArrayObjectType;

	typedef typename Superclass::DataObjectPointer                            DataObjectPointer;

	typedef ModelBase< InputImageType >                                       EnergyModelType;
	typedef typename EnergyModelType::Pointer                                 EnergyModelPointer;
	typedef typename EnergyModelType::ConstPointer                            EnergyModelConstPointer;

	typedef itk::ProcessObject ProcessObject;
	typedef itk::ThreadIdType ThreadIdType;
	typedef std::vector< MeasureArrayType >      ThreadMeasureArrayType;
	typedef std::vector< TotalVolumeContainer >  ThreadVolumeArrayType;

    void SetInput(const InputImageType *input) override {
    	this->SetNthInput(0, const_cast<InputImageType *>(input));
    }

    const InputImageType * GetInput() {
    	return static_cast<const InputImageType*>(this->ProcessObject::GetInput(0));
    }


    /** Set/Get priors
     *
     */
    void SetPriorsMap(const PriorsImageType *priors) {
    	this->SetNthInput(1, const_cast<PriorsImageType *>(priors));
    }

    const PriorsImageType * GetPriorsMap() {
    	return static_cast<const PriorsImageType*>(this->ProcessObject::GetInput(1));
    }

    /** Set/Get priors
     *
     */
    void SetMask(const MaskType *mask) {
    	this->SetNthInput(2, const_cast<MaskType *>(mask));
    }

    const MaskType * GetMask() {
    	return static_cast<const MaskType*>(this->ProcessObject::GetInput(2));
    }

    itkSetConstObjectMacro(Model, EnergyModelType);
    itkGetConstObjectMacro(Model, EnergyModelType);

	const MeasureArrayType GetEnergies() const { return this->GetEnergiesOutput()->Get();}
	MeasureArrayObjectType * GetEnergiesOutput();
	const MeasureArrayObjectType * GetEnergiesOutput() const;

protected:
	EnergyCalculatorFilter();
	virtual ~EnergyCalculatorFilter() {}
	void PrintSelf(std::ostream & os, itk::Indent indent) const override;

	void GenerateInputRequestedRegion() override;
	void BeforeThreadedGenerateData() override;
	void ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId) override;
	void AfterThreadedGenerateData() override;

	  /** Method that construct the outputs */
	typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
	using Superclass::MakeOutput;
	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) override;

private:
	EnergyCalculatorFilter(const Self &); //purposely not implemented
	void operator=(const Self &);         //purposely not implemented

	size_t m_NumberOfRegions;
	ThreadMeasureArrayType m_Energies;
	ThreadVolumeArrayType m_Volumes;
	PriorsPrecisionType m_PixelVolume;
	EnergyModelConstPointer m_Model;
}; // class EnergyCalculatorFilter


} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "EnergyCalculatorFilter.hxx"
#endif


#endif /* _ENERGYCALCULATORFILTER_H_ */
