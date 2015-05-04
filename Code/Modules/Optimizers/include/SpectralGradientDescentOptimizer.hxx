// --------------------------------------------------------------------------------------
// File:          SpectralGradientDescentOptimizer.hxx
// Date:          Jul 31, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of RegSeg
//
// RegSeg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RegSeg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RegSeg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef SPECTRALGRADIENTDESCENTOPTIMIZER_HXX_
#define SPECTRALGRADIENTDESCENTOPTIMIZER_HXX_

#include "SpectralGradientDescentOptimizer.h"

#include <itkImageAlgorithm.h>

using namespace std;

namespace rstk {

/**
 * Default constructor
 */
template< typename TFunctional >
SpectralGradientDescentOptimizer<TFunctional>::SpectralGradientDescentOptimizer() {}

template< typename TFunctional >
void SpectralGradientDescentOptimizer<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
}


template< typename TFunctional >
void SpectralGradientDescentOptimizer<TFunctional>::Iterate() {
	itkDebugMacro("Optimizer Iteration");
	this->ComputeUpdate(this->m_Coefficients, this->m_DerivativeCoefficients, this->m_NextCoefficients, true);
}

template< typename TFunctional >
void SpectralGradientDescentOptimizer<TFunctional>
::SetUpdate() {
	const typename CoefficientsImageType::PixelType* current[Dimension];
	for (size_t i = 0; i < Dimension; i++) {
		current[i] = this->m_NextCoefficients[i]->GetBufferPointer();
		itk::ImageAlgorithm::Copy< CoefficientsImageType, CoefficientsImageType > (
			this->m_NextCoefficients[i],
			this->m_Coefficients[i],
			this->m_NextCoefficients[i]->GetLargestPossibleRegion(),
			this->m_Coefficients[i]->GetLargestPossibleRegion()
		);
	}

	VectorType v;
	VectorType* buffer = this->m_CurrentCoefficients->GetBufferPointer();
	size_t nPix = this->m_NextCoefficients[0]->GetLargestPossibleRegion().GetNumberOfPixels();

	for(size_t i = 0; i < nPix; i++) {
		for(size_t d=0; d < Dimension; d++) {
			v[d] = *(current[d] + i);
		}
		*(buffer + i) = v;

	}
}

} // end namespace rstk

#endif /* SPECTRALGRADIENTDESCENTOPTIMIZER_HXX_ */
