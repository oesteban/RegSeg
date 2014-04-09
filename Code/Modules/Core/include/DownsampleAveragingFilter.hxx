/*
 * DownsampleAveragingFilter.hxx
 *
 *  Created on: Jun 26, 2013
 *      Author: oesteban
 */

#ifndef DOWNSAMPLEAVERAGINGFILTER_HXX_
#define DOWNSAMPLEAVERAGINGFILTER_HXX_

#include "DownsampleAveragingFilter.h"

#include <itkObjectFactory.h>
#include <itkIdentityTransform.h>
#include <itkProgressReporter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkSpecialCoordinatesImage.h>
#include <itkDefaultConvertPixelTraits.h>

namespace rstk
{
/**
 * Initialize new instance
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::DownsampleAveragingFilter()
{
  m_OutputOrigin.Fill(0.0);
  m_OutputSpacing.Fill(1.0);
  m_OutputDirection.SetIdentity();

  m_UseReferenceImage = false;

  m_Size.Fill(0);
  m_OutputStartIndex.Fill(0);

  m_DefaultPixelValue
    = itk::NumericTraits<PixelType>::ZeroValue( m_DefaultPixelValue );

  m_WindowN = 0;
}

/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "DefaultPixelValue: "
     << static_cast< typename itk::NumericTraits< PixelType >::PrintType >
  ( m_DefaultPixelValue )
     << std::endl;
  os << indent << "Size: " << m_Size << std::endl;
  os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
  os << indent << "UseReferenceImage: " << ( m_UseReferenceImage ? "On" : "Off" )
     << std::endl;
  return;
}

/**
 * Set the output image spacing.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::SetOutputSpacing(const double *spacing)
{
  SpacingType s;

  for( size_t i = 0; i < ImageDimension; i++ ) {
	  s = static_cast<SpacingValueType>(spacing[i]);
  }

  this->SetOutputSpacing(s);
}

/**
 * Set the output image origin.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::SetOutputOrigin(const double *origin)
{
  OriginPointType p(origin);

  this->SetOutputOrigin(p);
}

/** Helper method to set the output parameters based on this image */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::SetOutputParametersFromImage(const ImageBaseType *image)
{
  this->SetOutputOrigin ( image->GetOrigin() );
  this->SetOutputSpacing ( image->GetSpacing() );
  this->SetOutputDirection ( image->GetDirection() );
  this->SetOutputStartIndex ( image->GetLargestPossibleRegion().GetIndex() );
  this->SetSize ( image->GetLargestPossibleRegion().GetSize() );
}

/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::BeforeThreadedGenerateData()
{
  size_t tmpN = 1;

  SizeType inputSize = this->GetInput()->GetLargestPossibleRegion().GetSize();
  // Compute pixel ratios
  for( size_t i = 0; i<ImageDimension; i++) {
    m_WindowSize[i] = vcl_floor( inputSize[i] / m_Size[i] );
    if ( m_WindowSize[i] % 2 == 0 )
    	m_WindowSize[i]--;
    tmpN*=m_WindowSize[i];
  }
  m_WindowN = tmpN;

  unsigned int nComponents
    = itk::DefaultConvertPixelTraits<PixelType>::GetNumberOfComponents(
        m_DefaultPixelValue );

  if (nComponents == 0)
    {
    PixelComponentType zeroComponent
      = itk::NumericTraits<PixelComponentType>::ZeroValue( zeroComponent );
    nComponents = this->GetInput()->GetNumberOfComponentsPerPixel();
    itk::NumericTraits<PixelType>::SetLength(m_DefaultPixelValue, nComponents );
    for (unsigned int n=0; n<nComponents; n++)
      {
      PixelConvertType::SetNthComponent( n, m_DefaultPixelValue,
                                         zeroComponent );
      }
    }

}

/**
 * Set up state of filter after multi-threading.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::AfterThreadedGenerateData()
{

}

/**
 * ThreadedGenerateData
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  // Get the output pointers
  OutputImagePointer outputPtr = this->GetOutput();

  // Get this input pointers
  InputImageConstPointer inputPtr = this->GetInput();
  SizeType inputSize = inputPtr->GetLargestPossibleRegion().GetSize();

  // Create an iterator that will walk the output region for this thread.
  typedef itk::ImageRegionIteratorWithIndex< TOutputImage > OutputIterator;
  typedef itk::ImageRegionConstIterator< TInputImage > InputIterator;

  OutputIterator outIt(outputPtr, outputRegionForThread);

  // Define a few indices that will be used to translate from an input pixel
  // to an output pixel
  PointType outputPoint;         // Coordinates of current output pixel
  PointType inputPoint;          // Coordinates of current input pixel

  IndexType inputIndex;

  IndexType index;
  IndexType startInputIndex;
  size_t N;

  // Support for progress methods/callbacks
  itk::ProgressReporter progress( this,
                             threadId,
                             outputRegionForThread.GetNumberOfPixels() );

  // Cache information from the superclass
  PixelType defaultValue = this->GetDefaultPixelValue();

  InputImageRegionType window;
  SizeType windowSize;
  int offset[ImageDimension];
  for( size_t dim = 0; dim < ImageDimension; dim++) {
	  offset[dim] = - vcl_floor( 0.5*m_WindowSize[dim] );
  }

  while ( !outIt.IsAtEnd() )
    {
    // First get the position of the pixel in the output coordinate frame
    index = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint(index, outputPoint);

    // Compute corresponding input pixel continuous index, this index
    // will incremented in the scanline loop
    inputPtr->TransformPhysicalPointToIndex(outputPoint, inputIndex);

    PixelType pixval;
    double value = 0.0;

    windowSize = m_WindowSize;
    // Insert here averaging code
    for( size_t dim = 0; dim < ImageDimension; dim++) {
    	startInputIndex[dim] = inputIndex[dim] + offset[dim];

    	if ( startInputIndex[dim] < 0 ) {
    		windowSize[dim]+=startInputIndex[dim];
    		startInputIndex[dim] = 0;
    	}

    	if( (startInputIndex[dim] + windowSize[dim] ) > (inputSize[dim]-1) ) {
    		windowSize[dim] = inputSize[dim] - startInputIndex[dim] -1;
    	}

		if ( windowSize[dim]>m_WindowSize[dim] ) {
			windowSize[dim] = 0;
		}

    }

	window.SetSize( windowSize );
    window.SetIndex( startInputIndex );

    InputIterator inIt(inputPtr, window);

    N = 0;
    while( !inIt.IsAtEnd() ) {
    	value += inIt.Get();
    	N++;
    	++inIt;
    }

    pixval = (value!=0.0)?(value / N):0.0;
    outIt.Set(pixval);

    progress.CompletedPixel();
    ++outIt;
    } //while( !outIt.IsAtEnd() )

  return;
}

/**
 * Inform pipeline of necessary input image region
 *
 * Determining the actual input region is non-trivial, especially
 * when we cannot assume anything about the transform being used.
 * So we do the easy thing and request the entire input image.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
    return;
    }

  // get pointers to the input and output
  InputImagePointer inputPtr  =
    const_cast< TInputImage * >( this->GetInput() );

  // Request the entire input image
  inputPtr->SetRequestedRegionToLargestPossibleRegion();

  return;
}

/**
 * Set the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
const typename DownsampleAveragingFilter< TInputImage,
                                    TOutputImage,
                                    TPrecisionType >
::OutputImageType *
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::GetReferenceImage() const
{
  Self *                 surrogate = const_cast< Self * >( this );
  const OutputImageType *referenceImage =
    static_cast< const OutputImageType * >( surrogate->itk::ProcessObject::GetInput(1) );

  return referenceImage;
}

/**
 * Set the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::SetReferenceImage(const TOutputImage *image)
{
  itkDebugMacro("setting input ReferenceImage to " << image);
  if ( image != static_cast< const TOutputImage * >( this->itk::ProcessObject::GetInput(1) ) )
    {
    this->itk::ProcessObject::SetNthInput( 1, const_cast< TOutputImage * >( image ) );
    this->Modified();
    }
}

/**
 * Inform pipeline of required output region
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
void
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  OutputImagePointer outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }

  const OutputImageType *referenceImage = this->GetReferenceImage();

  // Set the size of the output region
  if ( m_UseReferenceImage && referenceImage )
    {
    outputPtr->SetLargestPossibleRegion(
      referenceImage->GetLargestPossibleRegion() );
    }
  else
    {
    typename TOutputImage::RegionType outputLargestPossibleRegion;
    outputLargestPossibleRegion.SetSize(m_Size);
    outputLargestPossibleRegion.SetIndex(m_OutputStartIndex);
    outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
    }

  // Set spacing and origin
  if ( m_UseReferenceImage && referenceImage )
    {
    outputPtr->SetSpacing( referenceImage->GetSpacing() );
    outputPtr->SetOrigin( referenceImage->GetOrigin() );
    outputPtr->SetDirection( referenceImage->GetDirection() );
    }
  else
    {
    outputPtr->SetSpacing(m_OutputSpacing);
    outputPtr->SetOrigin(m_OutputOrigin);
    outputPtr->SetDirection(m_OutputDirection);
    }

  return;
}

/**
 * Verify if any of the components has been modified.
 */
template< class TInputImage,
          class TOutputImage,
          class TPrecisionType >
itk::ModifiedTimeType
DownsampleAveragingFilter< TInputImage, TOutputImage, TPrecisionType >
::GetMTime(void) const
{
  itk::ModifiedTimeType latestTime = itk::Object::GetMTime();
  return latestTime;
}

} // end namespace rstk


#endif /* DOWNSAMPLEAVERAGINGFILTER_HXX_ */
