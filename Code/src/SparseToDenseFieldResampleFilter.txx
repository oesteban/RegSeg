/* --------------------------------------------------------------------------------------
 * File:    SparseToDenseFieldResampleFilter.txx
 * Date:    28/03/2012
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2012, Oscar Esteban - oesteban@die.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM) and
 Signal Processing Laboratory 5, EPFL (LTS5-EPFL).
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef SPARSETODENSEFIELDRESAMPLEFILTER_TXX_
#define SPARSETODENSEFIELDRESAMPLEFILTER_TXX_

#include "SparseToDenseFieldResampleFilter.h"
#include "IDWMultivariateInterpolator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace rstk {

template< class TInputMesh, class TOutputImage >
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >::
SparseToDenseFieldResampleFilter() {
	m_OutputSpacing.Fill(1.0);
	m_OutputOrigin.Fill(0.0);
	m_OutputDirection.SetIdentity();
	m_OutputSize.Fill(0);
	m_OutputStartIndex.Fill(0);

	m_Interpolator = IDWMultivariateInterpolator< InputMeshType, OutputImageType>::New();
	m_DefaultPixelValue = NumericTraits< ValueType >::ZeroValue( m_DefaultPixelValue );
};

/*
template< class TInputMesh, class TOutputImage >
void SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >::
GenerateData() {

};
*/
template< class TInputMesh, class TOutputImage >
void SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >::
PrintSelf( std::ostream &os, Indent indent) const {
	  //Superclass::PrintSelf(os, indent);

	  os << indent << "DefaultPixelValue: "
	     << static_cast< typename NumericTraits< OutputPixelType >::PrintType >( m_DefaultPixelValue )
	     << std::endl;
	  os << indent << "Size: " << m_OutputSize << std::endl;
	  os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
	  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
	  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
	  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
	  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
};


/**
 * Set the output image spacing.
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::SetOutputSpacing(const double *spacing)
{
  OutputSpacingType s(spacing);

  this->SetOutputSpacing(s);
}

/**
 * Set the output image origin.
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::SetOutputOrigin(const double *origin)
{
  OutputPointType p(origin);

  this->SetOutputOrigin(p);
}

/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::BeforeThreadedGenerateData()
{
  if ( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator not set");
    }

  // Connect input image to interpolator
  m_Interpolator->SetInputMesh( this->GetInput() );

  if (m_DefaultPixelValue.Size() == 0)
    {
    NumericTraits<OutputPixelType>::SetLength( m_DefaultPixelValue,
      this->GetOutput()->GetNumberOfComponentsPerPixel() );
    m_DefaultPixelValue.Fill(0);
    }
}

/**
 * Set up state of filter after multi-threading.
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::AfterThreadedGenerateData()
{
  // Disconnect input image from the interpolator
  m_Interpolator->SetInputImage(NULL);
}

/**
 * ThreadedGenerateData
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::ThreadedGenerateData( const OutputRegionType & outputRegionForThread,  ThreadIdType threadId )
{
  itkDebugMacro(<< "Currently executing");

  // Get the output pointers
  OutputImagePointer outputPtr = const_cast< TOutputImage * > (this->GetOutput());

  // Get ths input pointers
  InputMeshConstPointer inputPtr = this->GetInput();

  // Create an iterator that will walk the output region for this thread.
  typedef ImageRegionIteratorWithIndex< TOutputImage > OutputIterator;

  OutputIterator outIt(outputPtr, outputRegionForThread);

  // Define a few indices that will be used to translate from an input pixel
  // to an output pixel
  OutputPointType outputPoint;         // Coordinates of current output pixel
  OutputPointType inputPoint;          // Coordinates of current input pixel

  typedef ContinuousIndex< double, OutputImageDimension > ContinuousIndexType; // FIXME remove this double
  ContinuousIndexType inputIndex;

  // Doc says this only works for VectorImage, but Image implementation says otherwise...
  const unsigned int numberOfComponents = outputPtr->GetNumberOfComponentsPerPixel();

  // Support for progress methods/callbacks
  // FIXME
  //ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  typedef typename InterpolatorType::OutputType OutputType;

  // Walk the output region
  outIt.GoToBegin();


  OutputPixelType        pixval;
  NumericTraits< OutputPixelType >::SetLength( pixval, numberOfComponents );

  while ( !outIt.IsAtEnd() ) {
    outputPtr->TransformIndexToPhysicalPoint(outIt.GetIndex(), outputPoint);
    const OutputType value = m_Interpolator->Evaluate( outputPoint );

    for ( unsigned int i = 0; i < numberOfComponents; i++ ) {
      pixval[i] = static_cast< OutputVectorValueType >( value[i] );
    }

    outIt.Set(pixval);

    //progress.CompletedPixel();
    ++outIt;
  }
  return;
}

/**
 * Inform pipeline of necessary input image region
 *
 * Determining the actual input region is non-trivial, especially
 * when we cannot assume anything about the transform being used.
 * So we do the easy thing and request the entire input image.
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
    return;
    }

  // FIXME generateinputrequestedregion
/*
  // get pointers to the input and output
  OutputImagePointer inputPtr  =
    const_cast< TOutputImage * >( this->GetInput() );

  // Request the entire input image
  OutputRegionType inputRegion;
  inputRegion = inputPtr->GetLargestPossibleRegion();
  inputPtr->SetRequestedRegion(inputRegion);
*/
  return;
}

template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::CopyImageInformation( const OutputImageType* image ) {
	m_OutputDirection  = image->GetDirection();
	m_OutputOrigin     = image->GetOrigin();
	m_OutputSpacing    = image->GetSpacing();
	m_OutputStartIndex = image->GetRequestedRegion().GetIndex();
	m_OutputSize       = image->GetRequestedRegion().GetSize();
}

/**
 * Inform pipeline of required output region
 */
template< class TInputMesh, class TOutputImage >
void
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
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

  // Set the size of the output region
  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(m_OutputSize);
  outputLargestPossibleRegion.SetIndex(m_OutputStartIndex);
  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);

  // Set spacing and origin
  outputPtr->SetSpacing(m_OutputSpacing);
  outputPtr->SetOrigin(m_OutputOrigin);
  outputPtr->SetDirection(m_OutputDirection);

  return;
}

/**
 * Verify if any of the components has been modified.
 */
template< class TInputMesh, class TOutputImage >
unsigned long
SparseToDenseFieldResampleFilter< TInputMesh, TOutputImage >
::GetMTime(void) const
{
  unsigned long latestTime = Object::GetMTime();


  if ( m_Interpolator )
    {
    if ( latestTime < m_Interpolator->GetMTime() )
      {
      latestTime = m_Interpolator->GetMTime();
      }
    }

  return latestTime;
}

}

#endif /* SPARSETODENSEFIELDRESAMPLEFILTER_TXX_ */
