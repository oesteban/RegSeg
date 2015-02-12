// --------------------------------------------------------------------------------------
// File:          EnergyCalculatorFilter.h
// Date:          Dec 22, 2014
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

    void SetInput(const InputImageType *input) {
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
	void PrintSelf(std::ostream & os, itk::Indent indent) const;

	void GenerateInputRequestedRegion();
	void BeforeThreadedGenerateData();
	void ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId);
	void AfterThreadedGenerateData();

	  /** Method that construct the outputs */
	typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
	using Superclass::MakeOutput;
	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType);

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
