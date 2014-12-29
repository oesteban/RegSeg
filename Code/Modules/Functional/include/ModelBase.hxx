// --------------------------------------------------------------------------------------
// File:          ModelBase.hxx
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

#ifndef _MODELBASE_HXX_
#define _MODELBASE_HXX_

#include "ModelBase.h"

namespace rstk {
template< typename TInputVectorImage, typename TPriorsPrecisionType >
ModelBase< TInputVectorImage, TPriorsPrecisionType >
::ModelBase():
  Superclass(),
  m_NumberOfRegions(0),
  m_MaxEnergy(1.0e6) {
	m_MeasurementVectorSize = itk::NumericTraits<MeasurementVectorType>::GetLength(MeasurementVectorType());
	this->SetNumberOfRequiredInputs(3);
	this->SetNumberOfRequiredOutputs(1);
	this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
typename ModelBase< TInputVectorImage, TPriorsPrecisionType >::DataObjectPointer
ModelBase< TInputVectorImage, TPriorsPrecisionType >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed(idx) ) {
  return MembershipFunctionsObject::New().GetPointer();
}

} // namespace rstk
#endif /* _MODELBASE_HXX_ */
