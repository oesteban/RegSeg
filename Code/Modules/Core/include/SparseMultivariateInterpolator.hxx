// --------------------------------------------------------------------------
// File:             SparseMultivariateInterpolator.hxx
// Date:             30/10/2012
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
// This file is part of RegSeg
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

#ifndef SPARSEMULTIVARIATEINTERPOLATOR_HXX_
#define SPARSEMULTIVARIATEINTERPOLATOR_HXX_


#include "SparseMultivariateInterpolator.h"

namespace rstk {

/**
 * Constructor
 */
template<class TInputMesh, class TInputImage, class TCoordRep >
SparseMultivariateInterpolator< TInputMesh, TInputImage, TCoordRep >
::SparseMultivariateInterpolator():
 m_NumberOfImagePixels(0)
{
  m_Mesh = NULL;
}

/**
 * Standard "PrintSelf" method
 */
template<class TInputMesh, class TInputImage, class TCoordRep >
void
SparseMultivariateInterpolator< TInputMesh, TInputImage, TCoordRep >
::PrintSelf( std::ostream & os,  itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "InputMesh: " << m_Mesh.GetPointer() << std::endl;
}

/**
 * Initialize by setting the input image
 */
template<class TInputMesh, class TInputImage, class TCoordRep >
void
SparseMultivariateInterpolator< TInputMesh, TInputImage, TCoordRep >
::SetInputMesh( const InputMeshType *ptr)
{
  // set the input image
  m_Mesh = ptr;

  if ( ptr ) {
	  m_NumberOfMeshPoints = ptr->GetNumberOfPoints();
    //typename InputImageType::SizeType size = ptr->GetBufferedRegion().GetSize();
    //m_StartIndex = ptr->GetBufferedRegion().GetIndex();
    //
    //for ( unsigned int j = 0; j < ImageDimension; j++ )
    //  {
    //  m_EndIndex[j] = m_StartIndex[j] + static_cast< IndexValueType >( size[j] ) - 1;
    //  m_StartContinuousIndex[j] = static_cast< CoordRepType >( m_StartIndex[j] - 0.5 );
    //  m_EndContinuousIndex[j]   = static_cast< CoordRepType >( m_EndIndex[j] + 0.5 );
    //  }
  }
}

}



#endif /* SPARSEMULTIVARIATEINTERPOLATOR_HXX_ */
