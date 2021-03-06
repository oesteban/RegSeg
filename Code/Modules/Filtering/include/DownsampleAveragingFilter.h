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

#ifndef DOWNSAMPLEAVERAGINGFILTER_H_
#define DOWNSAMPLEAVERAGINGFILTER_H_


#include "itkFixedArray.h"
#include <itkTransform.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkExtrapolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSize.h>
#include <itkDefaultConvertPixelTraits.h>

namespace rstk
{
/** \class DownsampleAveragingFilter
 * \brief Shrink filter with averaging
 *
 * Output information (spacing, size and direction) for the output
 * image should be set. This information has the normal defaults of
 * unit spacing, zero origin and identity direction. Optionally, the
 * output information can be obtained from a reference image. If the
 * reference image is provided and UseReferenceImage is On, then the
 * spacing, origin and direction of the reference image will be used.
 *
 * Since this filter produces an image which is a different size than
 * its input, it needs to override several of the methods defined
 * in ProcessObject in order to properly manage the pipeline execution model.
 * In particular, this filter overrides
 * ProcessObject::GenerateInputRequestedRegion() and
 * ProcessObject::GenerateOutputInformation().
 *
 * This filter is implemented as a multithreaded filter.  It provides a
 * ThreadedGenerateData() method for its implementation.
 * \warning For multithreading, the TransformPoint method of the
 * user-designated coordinate transform must be threadsafe.
 *
 * \ingroup GeometricTransform
 * \ingroup ITKImageGrid
 *
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType = float>
class ITK_EXPORT DownsampleAveragingFilter :
  public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef DownsampleAveragingFilter                            Self;
  typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef itk::SmartPointer< Self >                            Pointer;
  typedef itk::SmartPointer< const Self >                      ConstPointer;

  typedef TInputImage                           InputImageType;
  typedef TOutputImage                          OutputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DownsampleAveragingFilter, ImageToImageFilter);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** base type for images of the current ImageDimension */
  typedef itk::ImageBase< itkGetStaticConstMacro(ImageDimension) > ImageBaseType;

  /** Image size typedef. */
  typedef itk::Size< itkGetStaticConstMacro(ImageDimension) > SizeType;

  /** Image index typedef. */
  typedef typename TOutputImage::IndexType IndexType;

  /** Image point typedef. */
  typedef typename TOutputImage::PointType            PointType;

  /** Image pixel value typedef. */
  typedef typename TOutputImage::PixelType PixelType;
  typedef typename TOutputImage::InternalPixelType OutputPixelValueType;
  typedef typename TInputImage::PixelType  InputPixelType;
  typedef typename TInputImage::InternalPixelType  InputPixelValueType;

  typedef itk::Image<OutputPixelValueType, ImageDimension > MaskImageType;

  typedef itk::DefaultConvertPixelTraits<PixelType> PixelConvertType;
  typedef itk::DefaultConvertPixelTraits<PixelType> InputPixelConvertType;

  typedef typename PixelConvertType::ComponentType PixelComponentType;

  /** Input pixel continuous index typdef */
  typedef itk::ContinuousIndex< TPrecisionType, ImageDimension >
  ContinuousInputIndexType;

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Image spacing,origin and direction typedef */
  typedef typename TOutputImage::SpacingType   SpacingType;
  typedef typename TOutputImage::SpacingValueType   SpacingValueType;
  typedef typename TOutputImage::PointType     OriginPointType;
  typedef typename TOutputImage::DirectionType DirectionType;

  /** Get/Set the size of the output image. */
  itkSetMacro(Size, SizeType);
  itkGetConstReferenceMacro(Size, SizeType);

  /** Get/Set the pixel value when a transformed pixel is outside of the
   * image.  The default default pixel value is 0. */
  itkSetMacro(DefaultPixelValue, PixelType);
  itkGetConstReferenceMacro(DefaultPixelValue, PixelType);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void SetOutputSpacing(const double *values);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, OriginPointType);
  virtual void SetOutputOrigin(const double *values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro(OutputOrigin, OriginPointType);

  /** Set the output direciton cosine matrix. */
  itkSetMacro(OutputDirection, DirectionType);
  itkGetConstReferenceMacro(OutputDirection, DirectionType);

  /** Helper method to set the output parameters based on this image */
  void SetOutputParametersFromImage(const ImageBaseType *image);

  /** Set the start index of the output largest possible region.
   * The default is an index of all zeros. */
  itkSetMacro(OutputStartIndex, IndexType);

  /** Get the start index of the output largest possible region. */
  itkGetConstReferenceMacro(OutputStartIndex, IndexType);

  /** Copy the output information from another Image.  By default,
   *  the information is specified with the SetOutputSpacing, Origin,
   *  and Direction methods. UseReferenceImage must be On and a
   *  Reference image must be present to override the defaul behavior.
   */
  void SetReferenceImage(const TOutputImage *image);
  const TOutputImage * GetReferenceImage() const;

  void SetMaskImage(const MaskImageType *image);
  const MaskImageType * GetMaskImage() const;

  itkSetMacro(UseReferenceImage, bool);
  itkBooleanMacro(UseReferenceImage);
  itkGetConstMacro(UseReferenceImage, bool);


  /** DownsampleAveragingFilter produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation() override;

  /** DownsampleAveragingFilter needs a different input requested region than
   * the output requested region.  As such, DownsampleAveragingFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() override;

  /** This method is used to set the state of the filter before
   * multi-threading. */
  virtual void BeforeThreadedGenerateData() override;

  /** This method is used to set the state of the filter after
   * multi-threading. */
  virtual void AfterThreadedGenerateData() override;

  /** Method Compute the Modified Time based on changed to the components. */
  itk::ModifiedTimeType GetMTime(void) const override;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputHasNumericTraitsCheck,
                   ( itk::Concept::HasNumericTraits< PixelComponentType > ) );
  /** End concept checking */
#endif

protected:
  DownsampleAveragingFilter();
  ~DownsampleAveragingFilter() {
  }
  void PrintSelf(std::ostream & os, itk::Indent indent) const override;

  /** Override VeriyInputInformation() since this filter's inputs do
   * not need to occoupy the same physical space.
   *
   * \sa ProcessObject::VerifyInputInformation
   */
  virtual void VerifyInputInformation() override {}

  /** DownsampleAveragingFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior
   * to calling ThreadedGenerateData().
   * ThreadedGenerateData can only write to the portion of the output image
   * specified by the parameter "outputRegionForThread"
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  virtual void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
		  itk::ThreadIdType threadId) override;

private:
  DownsampleAveragingFilter(const Self &); //purposely not implemented
  void operator=(const Self &);      //purposely not implemented

  SizeType        m_WindowSize;           // Size of the averaging window
  size_t          m_WindowN;
  PixelType       m_DefaultPixelValue;    // default pixel value
                                          // if the point is
                                          // outside the image
  PixelType       m_MaskedPixelValue;    // default pixel value
  SizeType        m_Size;                 // Size of the output image
  SpacingType     m_OutputSpacing;        // output image spacing
  OriginPointType m_OutputOrigin;         // output image origin
  DirectionType   m_OutputDirection;      // output image direction cosines
  IndexType       m_OutputStartIndex;     // output image start index
  size_t          m_NumberOfComponents;
  bool            m_UseReferenceImage;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "DownsampleAveragingFilter.hxx"
#endif

#endif /* DOWNSAMPLEAVERAGINGFILTER_H_ */
