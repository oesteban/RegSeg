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

#ifndef __VectorInterpolateImageFunction_h
#define __VectorInterpolateImageFunction_h

#include "itkImageFunction.h"
#include "itkFixedArray.h"

namespace rstk
{

/** \class VectorInterpolateImageFunction
 * \brief Base class for all vector image interpolaters.
 *
 * VectorInterpolateImageFunction is the base for all ImageFunctions that
 * interpolates image with vector pixel types. This function outputs
 * a return value of type Vector<double,Dimension>.
 *
 * This class is templated input image type and the coordinate
 * representation type.
 *
 * \warning This hierarchy of functions work only for images
 * with Vector-based pixel types. For scalar images use
 * InterpolateImageFunction.
 *
 * \sa InterpolateImageFunction
 * \ingroup ImageFunctions ImageInterpolators
 * \ingroup ITKImageFunction
 */
template< typename TInputImage, typename TCoordRep = double >
class VectorInterpolateImageFunction:
  public itk::ImageFunction<
    TInputImage,
    typename itk::NumericTraits< typename TInputImage::PixelType >::RealType,
    TCoordRep >
{
public:
  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef VectorInterpolateImageFunction Self;
  typedef itk::ImageFunction< TInputImage,
                         typename itk::NumericTraits< typename TInputImage::PixelType >::RealType,
                         TCoordRep >                          Superclass;

  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorInterpolateImageFunction, itk::ImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType           InputImageType;
  typedef typename InputImageType::PixelType            PixelType;
  typedef typename PixelType::ValueType                 ValueType;
  typedef typename itk::NumericTraits< ValueType >::RealType RealType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType      IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is RealType of TInputImage::PixelType. */
  typedef typename Superclass::OutputType OutputType;

  /** CoordRep typedef support. */
  typedef TCoordRep CoordRepType;

  /** Returns the interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType Evaluate(const PointType & point) const override
  {
    ContinuousIndexType index;

    this->GetInputImage()->TransformPhysicalPointToContinuousIndex(point, index);
    return ( this->EvaluateAtContinuousIndex(index) );
  }

  /** Interpolate the image at a continuous index position
   *
   * Returns the interpolated image intensity at a
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * Subclasses must override this method.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & index) const override = 0;

  /** Interpolate the image at an index position.
   * Simply returns the image value at the
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtIndex(const IndexType & index) const override
  {
    OutputType output;
    PixelType  input = this->GetInputImage()->GetPixel(index);

    for ( unsigned int k = 0;
              k < this->GetInputImage()->GetNumberOfComponentsPerPixel(); k++ )
      {
      output[k] = static_cast< double >( input[k] );
      }
    return ( output );
  }

protected:
  VectorInterpolateImageFunction() {}
  ~VectorInterpolateImageFunction() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const override
  { Superclass::PrintSelf(os, indent); }

private:
  VectorInterpolateImageFunction(const Self &); //purposely not implemented
  void operator=(const Self &);                 //purposely not implemented
};
} // end namespace itk

#endif
