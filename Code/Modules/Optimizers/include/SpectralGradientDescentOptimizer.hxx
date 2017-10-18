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
