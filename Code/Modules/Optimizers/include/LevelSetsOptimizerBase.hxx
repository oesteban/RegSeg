// --------------------------------------------------------------------------
// File:             LevelSetsOptimizerBase.hxx
// Date:             01/11/2012
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
// This file is part of ACWEReg
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

#ifndef LEVELSETSOPTIMIZERBASE_HXX_
#define LEVELSETSOPTIMIZERBASE_HXX_

#include "LevelSetsOptimizerBase.h"
#include "itkMultiThreader.h"



namespace rstk {

template< typename TLevelSetsFunction >
LevelSetsOptimizerBase<TLevelSetsFunction>::LevelSetsOptimizerBase() {
	this->m_LevelSets = NULL;
	this->m_CurrentLevelSetsValue = itk::NumericTraits<MeasureType>::infinity();
	this->m_ScalesAreIdentity = false;
}

template< typename TLevelSetsFunction >
void LevelSetsOptimizerBase<TLevelSetsFunction>::PrintSelf( std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);

	os << indent << "Number of threads: " << this->m_NumberOfThreads << std::endl;
	os << indent << "Number of scales:  " << this->m_Scales.Size() << std::endl;
//	if( this->m_Scales.Size() > 0 )	{
//		os << indent << "m_Scales: " << this->m_Scales << std::endl;
//	}
//	os << indent << "m_ScalesAreIdentity: " << this->GetScalesAreIdentity() << std::endl;
	os << indent << "LevelSets: " << std::endl;
	m_LevelSets->Print( os, indent.GetNextIndent() );
}

template< typename TLevelSetsFunction >
void LevelSetsOptimizerBase<TLevelSetsFunction>::Start() {
	/* Validate some settings */
	if( this->m_LevelSets.IsNull() )	{
		itkExceptionMacro("LevelSets object must be set.");
		return;
	}

	/* TODO Check A and B matrices */

}

template< typename TLevelSetsFunction >
const typename LevelSetsOptimizerBase<TLevelSetsFunction>::ParametersType &
LevelSetsOptimizerBase<TLevelSetsFunction>::GetCurrentPosition()
{
	if( this->m_LevelSets.IsNull() )
	{
		itkExceptionMacro("LevelSets has not been assigned. Cannot get parameters.");
	}
	return this->m_LevelSets->GetParameters();
}

//-------------------------------------------------------------------
template< typename TLevelSetsFunction >
const typename LevelSetsOptimizerBase<TLevelSetsFunction>::MeasureType &
LevelSetsOptimizerBase<TLevelSetsFunction>::GetValue(){
	return this->GetCurrentLevelSetsValue();
}

}


#endif /* LEVELSETSOPTIMIZERBASE_HXX_ */
