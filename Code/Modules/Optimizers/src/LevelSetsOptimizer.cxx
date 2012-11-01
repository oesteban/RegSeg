// --------------------------------------------------------------------------
// File:             LevelSetsOptimizer.txx
// Date:             15/10/2012
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
// This file is part of ACWE-Registration
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

#include "LevelSetsOptimizer.h"

namespace rstk {

LevelSetsOptimizer::LevelSetsOptimizer() {
	/* Initialize state tracking variables */
	this->m_NumberOfIterations = 100;
	this->m_CurrentIteration   = 0;
	this->m_StopCondition      = MAXIMUM_NUMBER_OF_ITERATIONS;
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";

}

LevelSetsOptimizer::~LevelSetsOptimizer() { }

void LevelSetsOptimizer::PrintSelf( std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
	os << indent << "Number of iterations: " << this->m_NumberOfIterations << std::endl;
	os << indent << "Current iteration: " << this->m_CurrentIteration << std::endl;
	os << indent << "Stop condition:" << this->m_StopCondition << std::endl;
	os << indent << "Stop condition description: " << this->m_StopConditionDescription.str() << std::endl;
}

const LevelSetsOptimizer::StopConditionReturnStringType
LevelSetsOptimizer::GetStopConditionDescription() const {
  return this->m_StopConditionDescription.str();
}

void LevelSetsOptimizer::Stop(void) {
  itkDebugMacro( "StopOptimization called with a description - "
    << this->GetStopConditionDescription() );
  this->m_Stop = true;
  this->InvokeEvent( itk::EndEvent() );
}

} // End of namespace rstk


