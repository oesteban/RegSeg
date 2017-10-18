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

#ifndef _ENERGYCALCULATORFILTER_HXX_
#define _ENERGYCALCULATORFILTER_HXX_

#include "EnergyCalculatorFilter.h"
#include <itkImageBase.h>
#include <itkProgressReporter.h>

namespace rstk {

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::EnergyCalculatorFilter():
 Superclass(),
 m_PixelVolume(1.0) {
	this->SetNumberOfRequiredInputs(3);
	this->SetNumberOfRequiredOutputs(1);
	this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );
};

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
typename EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >::DataObjectPointer
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) ) {
  return MeasureArrayObjectType::New().GetPointer();
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::BeforeThreadedGenerateData() {
	// find the actual number of threads
	long nbOfThreads = this->GetNumberOfThreads();
	if ( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 ) {
		nbOfThreads = vnl_math_min( this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
	}
	// number of threads can be constrained by the region size, so call the
	// SplitRequestedRegion
	// to get the real number of threads which will be used
	RegionType splitRegion;  // dummy region - just to call the following method
	nbOfThreads = this->SplitRequestedRegion(0, nbOfThreads, splitRegion);

	this->m_NumberOfRegions = this->GetPriorsMap()->GetNumberOfComponentsPerPixel();
	this->m_Energies.resize(nbOfThreads);
	this->m_Volumes.resize(nbOfThreads);

	this->m_PixelVolume = 1.0;
	SpacingType s = this->GetInput()->GetSpacing();
	for(size_t i = 0; i<Dimension; i++) this->m_PixelVolume*= s[i];
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId) {
	long nbOfPixels = inputRegionForThread.GetNumberOfPixels();
	itk::ProgressReporter progress( this, threadId, nbOfPixels );

	itk::ImageRegionConstIterator< TInputVectorImage > inputIt( this->GetInput(), inputRegionForThread );
	itk::ImageRegionConstIterator< PriorsImageType >   priorIt( this->GetPriorsMap(), inputRegionForThread );
	itk::ImageRegionConstIterator< MaskType >          maskIt ( this->GetMask(), inputRegionForThread );
 	inputIt.GoToBegin();

 	EnergyModelConstPointer model = this->GetModel();

 	MeasureArrayType energies;
 	energies.SetSize(this->m_NumberOfRegions);
 	energies.Fill(0.0);

 	TotalVolumeContainer volumes;
 	volumes.SetSize(this->m_NumberOfRegions);
 	volumes.Fill(0.0);

 	size_t lastroi = this->m_NumberOfRegions - 1;
 	PriorsPixelType w;
 	double bgw;
 	PriorsPrecisionType vol;
 	PixelType val;
 	MeasureType e;
	while ( !inputIt.IsAtEnd() ) {
		bgw = maskIt.Get();
		val = inputIt.Get();
		w = priorIt.Get();

		for(size_t roi = 0; roi < m_NumberOfRegions; roi++ ) {
			if( w[roi] < 1.0e-8 )
				continue;

			e = model->Evaluate(val, roi);
			vol = w[roi] * this->m_PixelVolume;
			volumes[roi]+= vol;
			energies[roi]+= vol * e;
		}

		++inputIt;
	    ++priorIt;
	    progress.CompletedPixel();
	}

	this->m_Energies[threadId] = energies;
	this->m_Volumes[threadId] = volumes;
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::AfterThreadedGenerateData() {
 	MeasureArrayType energies;
 	energies.SetSize(this->m_NumberOfRegions);
 	energies.Fill(0.0);

 	TotalVolumeContainer volumes;
 	volumes.SetSize(this->m_NumberOfRegions);
 	volumes.Fill(0.0);

	EnergyModelConstPointer model = this->GetModel();
	for(size_t roi = 0; roi < m_NumberOfRegions; roi++ ) {

		for(size_t th = 0; th < this->m_Volumes.size(); th++) {
			volumes[roi]+= this->m_Volumes[th][roi];
			energies[roi]+= this->m_Energies[th][roi];
		}
		energies[roi]+= volumes[roi] * model->GetRegionOffsetContainer()[roi];
	}

	this->GetEnergiesOutput()->Set(energies);
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::GenerateInputRequestedRegion()
 {

 }


template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::PrintSelf(std::ostream & os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
typename EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >::MeasureArrayObjectType *
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::GetEnergiesOutput()
{
	return static_cast<MeasureArrayObjectType *>(this->ProcessObject::GetOutput(0));
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
const typename EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >::MeasureArrayObjectType *
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::GetEnergiesOutput() const
{
	return static_cast<const MeasureArrayObjectType *>(this->ProcessObject::GetOutput(0));
}

//template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
//typename EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >::MeasureArrayType
//EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
//::GetOutput()
//{
//  const MeasureArrayType *output =
//    static_cast< const MeasureArrayType * >( this->ProcessObject::GetOutput(0));
//  return output;
//}

} // namespace rstk
#endif /* _ENERGYCALCULATORFILTER_HXX_ */
