// --------------------------------------------------------------------------
// File:             EnergyOptimizerBase.cxx
// Date:             17/10/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWERegistration-Debug@Debug
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "EnergyOptimizerBase.h"
#include "itkMultiThreader.h"



namespace rstk {

EnergyOptimizerBase::EnergyOptimizerBase() {
	this->m_Energy = NULL;
	this->m_CurrentEnergyValue = itk::NumericTraits<MeasureType>::infinity();
	this->m_ScalesAreIdentity = false;
}

EnergyOptimizerBase::~EnergyOptimizerBase() {

}


void EnergyOptimizerBase::PrintSelf( std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);

	os << indent << "Number of threads: " << this->m_NumberOfThreads << std::endl;
	os << indent << "Number of scales:  " << this->m_Scales.Size() << std::endl;
//	if( this->m_Scales.Size() > 0 )	{
//		os << indent << "m_Scales: " << this->m_Scales << std::endl;
//	}
//	os << indent << "m_ScalesAreIdentity: " << this->GetScalesAreIdentity() << std::endl;
	os << indent << "Energy: " << std::endl;
	m_Energy->Print( os, indent.GetNextIndent() );
}


void EnergyOptimizerBase::Start() {
	/* Validate some settings */
	if( this->m_Energy.IsNull() )	{
		itkExceptionMacro("Energy object must be set.");
		return;
	}

	/* Verify m_Scales. If m_Scales hasn't been set, initialize to all 1's. */
	typedef ScalesType::ValueType     ValueType;
	if( this->m_Scales.Size() > 0 ) {
		if( this->m_Scales.Size() != this->m_Energy->GetNumberOfLocalParameters() ) {
			itkExceptionMacro("Size of scales (" << this->m_Scales.Size()
					<< ") must equal number of local parameters (" <<
					this->m_Energy->GetNumberOfLocalParameters() << ").");
		}
		/* Check that all values in m_Scales are > machine epsilon, to avoid
		 * division by zero/epsilon.
		 * Also check if scales are identity. */
		typedef ScalesType::size_type     SizeType;
		this->m_ScalesAreIdentity = true;
		for( SizeType i=0; i < this->m_Scales.Size(); i++ ) {
			if( this->m_Scales[i] <= itk::NumericTraits<ValueType>::epsilon() ) {
				itkExceptionMacro("m_Scales values must be > epsilon.");
			}
			/* Check if the scales are identity. Consider to be identity if
			 * within a tolerance, to allow for automatically estimated scales
			 * that may not be exactly 1.0 when in priciniple they should be. */
			ValueType difference =
					vcl_fabs( itk::NumericTraits<ValueType>::OneValue() - this->m_Scales[i] );
			ValueType tolerance = static_cast<ValueType>( 0.01 );
			if( difference > tolerance  ) {
				this->m_ScalesAreIdentity = false;
				break;
			}
		}
	}
	else {
		m_Scales.SetSize( this->m_Energy->GetNumberOfLocalParameters() );
		m_Scales.Fill( itk::NumericTraits<ValueType>::OneValue() );
		this->m_ScalesAreIdentity = true;
	}
}

const EnergyOptimizerBase::ParametersType &
EnergyOptimizerBase::GetCurrentPosition()
{
	if( this->m_Energy.IsNull() )
	{
		itkExceptionMacro("Energy has not been assigned. Cannot get parameters.");
	}
	return this->m_Energy->GetParameters();
}

//-------------------------------------------------------------------
const EnergyOptimizerBase::MeasureType &
EnergyOptimizerBase::GetValue(){
	return this->GetCurrentEnergyValue();
}

}
