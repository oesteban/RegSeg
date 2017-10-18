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

#ifndef __VectorLinearInterpolateImageFunction_hxx
#define __VectorLinearInterpolateImageFunction_hxx

#include "VectorLinearInterpolateImageFunction.h"

#include "vnl/vnl_math.h"

namespace rstk
{
/**
 * Define the number of neighbors
 */
template< typename TInputImage, typename TCoordRep >
const unsigned long
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;

/**
 * Constructor
 */
template< typename TInputImage, typename TCoordRep >
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::VectorLinearInterpolateImageFunction()
{}

/**
 * PrintSelf
 */
template< typename TInputImage, typename TCoordRep >
void
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

/**
 * Evaluate at image index position
 */
template< typename TInputImage, typename TCoordRep >
typename VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
VectorLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType & index) const
{
  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  IndexType baseIndex;
  InternalComputationType    distance[ImageDimension];
  const TInputImage * const inputImgPtr=this->GetInputImage();
  size_t ncomps = inputImgPtr->GetNumberOfComponentsPerPixel();

  for (unsigned int dim = 0; dim < ImageDimension; ++dim )
    {
    baseIndex[dim] = itk::Math::Floor< itk::IndexValueType >(index[dim]);
    distance[dim] = index[dim] - static_cast< InternalComputationType >( baseIndex[dim] );
    }

  /**
   * Interpolated value is the weighted sum of each of the surrounding
   * neighbors. The weight for each neighbor is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  OutputType output;
  itk::NumericTraits<OutputType>::SetLength(output, ncomps);
  output.Fill(0.0);

  typedef typename itk::NumericTraits< PixelType >::ScalarRealType ScalarRealType;
  ScalarRealType totalOverlap = itk::NumericTraits< ScalarRealType >::Zero;

  for ( unsigned int counter = 0; counter < m_Neighbors; ++counter )
    {
    InternalComputationType overlap = 1.0;    // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour

    IndexType    neighIndex;
    // get neighbor index and overlap fraction
    for ( unsigned int dim = 0; dim < ImageDimension; ++dim )
      {
      if ( upper & 1 )
        {
        neighIndex[dim] = baseIndex[dim] + 1;
        // Take care of the case where the pixel is just
        // in the outer upper boundary of the image grid.
        if ( neighIndex[dim] > this->m_EndIndex[dim] )
          {
          neighIndex[dim] = this->m_EndIndex[dim];
          }
        overlap *= distance[dim];
        }
      else
        {
        neighIndex[dim] = baseIndex[dim];
        // Take care of the case where the pixel is just
        // in the outer lower boundary of the image grid.
        if ( neighIndex[dim] < this->m_StartIndex[dim] )
          {
          neighIndex[dim] = this->m_StartIndex[dim];
          }
        overlap *= 1.0 - distance[dim];
        }
      upper >>= 1;
      }

    // get neighbor value only if overlap is not zero
    if ( overlap )
      {
      const PixelType & input = inputImgPtr->GetPixel(neighIndex);
      for ( unsigned int k = 0; k < ncomps; ++k )
        {
        output[k] += overlap * static_cast< InternalComputationType >( input[k] );
        }
      totalOverlap += overlap;
      }

    if ( totalOverlap == 1.0 )
      {
      // finished
      break;
      }
    }

  return ( output );
}
} // end namespace itk

#endif
