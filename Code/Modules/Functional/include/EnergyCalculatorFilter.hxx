// --------------------------------------------------------------------------------------
// File:          EnergyCalculatorFilter.hxx
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

#ifndef _ENERGYCALCULATORFILTER_HXX_
#define _ENERGYCALCULATORFILTER_HXX_

#include "EnergyCalculatorFilter.h"
#include <itkImageBase.h>

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
	this->m_Energies.SetSize(this->m_NumberOfRegions);
	this->m_Energies.Fill(0.0);
	this->m_Volumes.SetSize(this->m_NumberOfRegions);
	this->m_Volumes.Fill(0.0);

	this->m_PixelVolume = 1.0;
	SpacingType s = this->GetInput()->GetSpacing();
	for(size_t i = 0; i<Dimension; i++) this->m_PixelVolume*= s[i];
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::ThreadedGenerateData(const RegionType & inputRegionForThread, ThreadIdType threadId) {
	long nbOfPixels = inputRegionForThread.GetNumberOfPixels();
	ProgressReporter progress( this, threadId, nbOfPixels );

	ImageRegionConstIterator< TInputVectorImage > inputIt( this->GetInput(), inputRegionForThread );
	ImageRegionConstIterator< PriorsImageType >   priorIt( this->GetPriorsMap(), inputRegionForThread );
	ImageRegionConstIterator< MaskType >          maskIt ( this->GetMask(), inputRegionForThread );
 	inputIt.GoToBegin();

 	EnergyModelConstPointer model = this->GetModel();

 	PriorsPixelType w;
 	PriorsPrecisionType m, vol;
 	PixelType v;
	while ( !inputIt.IsAtEnd() ) {
		m = maskIt.Get();

		if (m > 1.0e-5) {
			v = inputIt.Get();
			w = priorIt.Get();

			for(size_t roi = 0; roi < m_NumberOfRegions; roi++ ) {
				vol = w[roi] * this->m_PixelVolume;
				this->m_Volumes[roi]+= vol;
				this->m_Energies[roi]+= vol * model->Evaluate(v, roi);
			}
		}

	    ++inputIt;
	    ++priorIt;
	    progress.CompletedPixel();
	}
}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::AfterThreadedGenerateData() {
	EnergyModelConstPointer model = this->GetModel();


	std::cout << this->m_Energies;
	for(size_t roi = 0; roi < m_NumberOfRegions; roi++ ) {
		this->m_Energies[roi]+= this->m_Volumes[roi] * model->GetRegionOffsetContainer()[roi];
	}

	std::cout << " w/offsets=" << this->m_Energies << std::endl;

}

template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::GenerateInputRequestedRegion()
 {
	this->m_Energies.SetSize( itk::NumericTraits<PixelType>::GetLength( PixelType() ) );
 }


template < typename TInputVectorImage, typename TMeasureType, typename TPriorsPrecisionType >
void
EnergyCalculatorFilter< TInputVectorImage, TMeasureType, TPriorsPrecisionType >
::PrintSelf(std::ostream & os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
}

}

#endif /* _ENERGYCALCULATORFILTER_HXX_ */
