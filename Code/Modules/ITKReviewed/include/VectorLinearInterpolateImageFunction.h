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

#ifndef __VectorLinearInterpolateImageFunction_h
#define __VectorLinearInterpolateImageFunction_h

#include "VectorInterpolateImageFunction.h"

namespace rstk
{
/**
 * \class VectorLinearInterpolateImageFunction
 * \brief Linearly interpolate a vector image at specified positions.
 *
 * VectorLinearInterpolateImageFunction linearly interpolates a vector
 * image intensity non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type.
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for Vector images. For
 * scalar images use LinearInterpolateImageFunction.
 *
 * \ingroup ImageFunctions ImageInterpolators
 *
 * \ingroup ITKImageFunction
 */
template< typename TInputImage, typename TCoordRep = double >
class VectorLinearInterpolateImageFunction:
  public VectorInterpolateImageFunction< TInputImage, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef VectorLinearInterpolateImageFunction                     Self;
  typedef VectorInterpolateImageFunction< TInputImage, TCoordRep > Superclass;
  typedef itk::SmartPointer< Self >                                     Pointer;
  typedef itk::SmartPointer< const Self >                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorLinearInterpolateImageFunction,
		       VectorInterpolateImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::PixelType      PixelType;
  typedef typename Superclass::ValueType      ValueType;
  typedef typename Superclass::RealType       RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType      IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
  typedef typename ContinuousIndexType::ValueType  InternalComputationType;

  /** Output type is Vector<double,Dimension> */
  typedef typename Superclass::OutputType OutputType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & index) const override;

protected:
  VectorLinearInterpolateImageFunction();
  ~VectorLinearInterpolateImageFunction(){}
  void PrintSelf(std::ostream & os, itk::Indent indent) const override;

private:
  VectorLinearInterpolateImageFunction(const Self &); //purposely not
                                                      // implemented
  void operator=(const Self &);                       //purposely not
                                                      // implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long m_Neighbors;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "VectorLinearInterpolateImageFunction.hxx"
#endif

#endif
