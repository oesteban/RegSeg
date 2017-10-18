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

#ifndef __itkInternalOrientationFilter_h
#define __itkInternalOrientationFilter_h

#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkSpatialOrientationAdapter.h"
#include <map>
#include <string>

namespace rstk
{
/** \class InternalOrientationFilter
 * \brief Permute axes and then flip images as needed to obtain
 *  agreement in coordinateOrientation codes.
 *
 * This class satisfies a common requirement in medical imaging, which
 * is to properly orient a 3 dimensional image with respect to anatomical
 * features.  Due to the wide variety of hardware used to generate 3D images
 * of human anatomy, and the even wider variety of image processing software,
 * it is often necessary to re-orient image volume data.
 *
 * InternalOrientationFilter depends on a set of constants that describe all possible
 * labeled according to the following scheme:
 * Directions are labeled in terms of following pairs:
 *   - Left and Right (Subject's left and right)
 *   - Anterior and Posterior (Subject's front and back)
 *   - Inferior and Superior (Subject's bottom and top, i.e. feet and head)
 *
 * The initials of these directions are used in a 3 letter code in the
 * enumerated type itk::SpatialOrientation::ValidCoordinateOrientationFlags.
 * The initials are given fastest moving index first, second fastest second,
 * third fastest third.
 * Examples:
 *  - ITK_COORDINATE_ORIENTATION_RIP
 *    -# Right to Left varies fastest (0th pixel on Subject's right)
 *    -# Inferior to Superior varies second fastest
 *    -# Posterior to Anterior varies slowest.
 *  - ITK_COORDINATE_ORIENTATION_LSA
 *    -# Left to Right varies fastest (0th pixel on Subject's left)
 *    -# Superior to Inferior varies second fastest
 *    -# Anterior to Posterior varies slower
 *
 * In order to use this filter, you need to supply an input
 * image, the current orientation of the input image (set with
 * SetGivenCoordinateOrientation) and the desired orientation.
 * (set with SetDesiredCoordinateOrientation).
 * You may explicitly set the DesiredOrientation with
 * SetDesiredCoordinateOrientation (if UseImageDirection is "off") or
 * you can use the image's direction cosines to set the
 * DesiredOrientation (if UseImageDirection is "on").
 * When reading image files that define the coordinate orientation
 * of the image, the current orientation is stored in the MetadataDictionary
 * for the itk::Image object and the Image.Direction direction cosine
 * matrix created from the file.
 *
 * As an example, if you wished to keep all images within your program in the
 * orientation corresponding to the Analyze file format's 'CORONAL' orientation
 * you could do something like the following
 *
 * \code
 * // DEPRECATED -- using metadata for orientation is no longer supported
 * //
 * #include "itkAnalyzeImageIO.h"
 * #include "itkMetaDataObject.h"
 * #include "itkImage.h"
 * #include "itkInternalOrientationFilter.h"
 * typedef itk::Image<unsigned char,3> ImageType;
 * typedef itk::ImageFileReader< TstImageType > ImageReaderType;
 * ImageType::Pointer ReadAnalyzeFile(const char *path)
 * {
 *   itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
 *   ImageReaderType::Pointer fileReader = ImageReaderType::New();
 *   fileReader->SetImageIO(io);
 *   fileReader->SetFileName(path);
 *   fileReader->Update();
 *   ImageType::Pointer rval = fileReader->GetOutput();
 *
 * Deprecated -- use direction cosines
 *  itk::SpatialOrientation::ValidCoordinateOrientationFlags fileOrientation;
 *  itk::ExposeMetaData<itk::SpatialOrientation::ValidCoordinateOrientationFlags>
 *    (rval->GetMetaDataDictionary(),itk::ITK_CoordinateOrientation,fileOrientation);
 *   itk::InternalOrientationFilter<ImageType,ImageType>::Pointer orienter =
 *     itk::InternalOrientationFilter<ImageType,ImageType>::New();
 *   orienter->SetGivenCoordinateOrientation(fileOrientation); // deprecated
 *
 *   orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP);
 *   orienter->SetInput(rval);
 *   orienter->Update();
 *   rval = orienter->GetOutput();
 *   return rval;
 * }
 * \endcode
 *
 * Or, using the direction cosines of the image,
 * \code
 * #include "itkAnalyzeImageIO.h"
 * #include "itkImage.h"
 * #include "itkInternalOrientationFilter.h"
 * typedef itk::Image<unsigned char,3> ImageType;
 * typedef itk::ImageFileReader< TstImageType > ImageReaderType;
 * ImageType::Pointer ReadAnalyzeFile(const char *path)
 * {
 *   itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
 *   ImageReaderType::Pointer fileReader = ImageReaderType::New();
 *   fileReader->SetImageIO(io);
 *   fileReader->SetFileName(path);
 *   fileReader->Update();
 *   ImageType::Pointer rval = fileReader->GetOutput();
 *
 *   itk::InternalOrientationFilter<ImageType,ImageType>::Pointer orienter =
 *     itk::InternalOrientationFilter<ImageType,ImageType>::New();
 *   orienter->UseImageDirectionOn();
 *   orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP);
 *   orienter->SetInput(rval);
 *   orienter->Update();
 *   rval = orienter->GetOutput();
 *   return rval;
 * }
 * \endcode
 * \ingroup ITKImageGrid
 */
template< typename TInputImage, typename TOutputImage >
class InternalOrientationFilter:
  public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef InternalOrientationFilter                               Self;
  typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef itk::SmartPointer< Self >                            Pointer;
  typedef itk::SmartPointer< const Self >                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename InputImageType::Pointer       InputImagePointer;
  typedef typename InputImageType::ConstPointer  InputImageConstPointer;
  typedef typename InputImageType::RegionType    InputImageRegionType;
  typedef typename InputImageType::PixelType     InputImagePixelType;
  typedef typename InputImageType::PointType     InputImagePointType;
  typedef typename InputImageType::DirectionType DirectionType;
  typedef TOutputImage                           OutputImageType;
  typedef typename OutputImageType::Pointer      OutputImagePointer;
  typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
  typedef typename OutputImageType::RegionType   OutputImageRegionType;
  typedef typename OutputImageType::PixelType    OutputImagePixelType;
  typedef itk::SpatialOrientation::ValidCoordinateOrientationFlags
  CoordinateOrientationCode;
  /** Axes permuter type. */
  typedef itk::PermuteAxesImageFilter< TInputImage >        PermuterType;
  typedef typename PermuterType::PermuteOrderArrayType PermuteOrderArrayType;

  /** Axes flipper type. */
  typedef itk::FlipImageFilter< TInputImage >          FlipperType;
  typedef typename FlipperType::FlipAxesArrayType FlipAxesArrayType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(InternalOrientationFilter, itk::ImageToImageFilter);

  /** Get axes permute order. */
  itkGetConstReferenceMacro(PermuteOrder, PermuteOrderArrayType);

  /** Get flip axes. */
  itkGetConstReferenceMacro(FlipAxes, FlipAxesArrayType);


  /** InternalOrientationFilter produces an image which is a different
   * dimensionality than its input image, in general. As such,
   * InternalOrientationFilter needs to provide an implementation for
   * GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below.
   * \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation();

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputConvertibleToOutput,
                   ( itk::Concept::Convertible< InputImagePixelType, OutputImagePixelType > ) );
  itkConceptMacro( SameDimension,
                   ( itk::Concept::SameDimension< itkGetStaticConstMacro(InputImageDimension),
                                             itkGetStaticConstMacro(OutputImageDimension) > ) );
  itkConceptMacro( DimensionShouldBe3,
                   ( itk::Concept::SameDimension< itkGetStaticConstMacro(InputImageDimension), 3 > ) );
  // End concept checking
#endif

protected:
  InternalOrientationFilter();
  ~InternalOrientationFilter() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

  /** InternalOrientationFilter needs the entire input be
   * available. Thus, it needs to provide an implementation of
   * GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** InternalOrientationFilter will produce the entire output. */
  void EnlargeOutputRequestedRegion( itk::DataObject *itkNotUsed(output) );

  /*** Member functions used by GenerateData: */
  void DeterminePermutationsAndFlips(const itk::SpatialOrientation::ValidCoordinateOrientationFlags fixed_orient,
                                     const itk::SpatialOrientation::ValidCoordinateOrientationFlags moving_orient);

  bool NeedToPermute();

  bool NeedToFlip();

  /** Single-threaded version of GenerateData.  This filter delegates
   * to PermuteAxesImageFilter and FlipImageFilter. */
  void GenerateData();

private:
  InternalOrientationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);    //purposely not implemented

  std::string GetMajorAxisFromPatientRelativeDirectionCosine(double x, double y, double z);

  CoordinateOrientationCode m_GivenCoordinateOrientation;
  CoordinateOrientationCode m_DesiredCoordinateOrientation;
  bool                      m_UseImageDirection;

  PermuteOrderArrayType m_PermuteOrder;
  FlipAxesArrayType     m_FlipAxes;

  std::map< std::string, CoordinateOrientationCode > m_StringToCode;
  std::map< CoordinateOrientationCode, std::string > m_CodeToString;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "InternalOrientationFilter.hxx"
#endif

#endif
