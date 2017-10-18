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

#ifndef _MODELBASE_HXX_
#define _MODELBASE_HXX_

#include "ModelBase.h"

namespace rstk {
template< typename TInputVectorImage, typename TPriorsPrecisionType >
ModelBase< TInputVectorImage, TPriorsPrecisionType >
::ModelBase():
  Superclass(),
  m_NumberOfRegions(0),
  m_MaxEnergy(0.0),
  m_NumberOfSpecialRegions(1) {
	m_MeasurementVectorSize = itk::NumericTraits<MeasurementVectorType>::GetLength(MeasurementVectorType());
	this->SetNumberOfRequiredInputs(3);
	this->SetNumberOfRequiredOutputs(1);
	this->SetNthOutput( 0, this->MakeOutput(0) );

}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
typename ModelBase< TInputVectorImage, TPriorsPrecisionType >::DataObjectPointer
ModelBase< TInputVectorImage, TPriorsPrecisionType >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) ) {
  return MembershipFunctionsObject::New().GetPointer();
}

} // namespace rstk
#endif /* _MODELBASE_HXX_ */
