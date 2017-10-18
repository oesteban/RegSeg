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
