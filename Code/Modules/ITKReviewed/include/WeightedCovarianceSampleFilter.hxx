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

#ifndef __WeightedCovarianceSampleFilter_hxx
#define __WeightedCovarianceSampleFilter_hxx

#include "WeightedCovarianceSampleFilter.h"
#include <itkWeightedMeanSampleFilter.h>

#include <boost/math/special_functions/digamma.hpp>

namespace bm = boost::math;

namespace itk
{
template< typename TSample >
WeightedCovarianceSampleFilter< TSample >
::WeightedCovarianceSampleFilter():
 m_RemoveOutliers(true),
 Superclass()
{
  this->ProcessObject::SetNthInput(1, ITK_NULLPTR);

  this->ProcessObject::SetNumberOfRequiredOutputs(4);
  this->ProcessObject::SetNthOutput( 2, this->MakeOutput(2) );
  this->ProcessObject::SetNthOutput( 3, this->MakeOutput(3) );
}

template< typename TSample >
WeightedCovarianceSampleFilter< TSample >
::~WeightedCovarianceSampleFilter()
{}

template< typename TSample >
void
WeightedCovarianceSampleFilter< TSample >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  // m_Weights
  os << indent << "Weights: " << this->GetWeightsInput() << std::endl;
}

template< typename TSample >
inline void
WeightedCovarianceSampleFilter< TSample >
::GenerateData()
{
  // if weight array is specified use it to compute the covariance
  const InputWeightArrayObjectType *weightArrayObject =
    this->GetWeightsInput();

  if ( weightArrayObject != ITK_NULLPTR )
    {
    this->ComputeCovarianceMatrixWithWeights();
    return;
    }

  // Otherwise compute the regular covariance matrix ( without weight
  // coefficients)
  Superclass::GenerateData();
}


template< typename TSample >
typename WeightedCovarianceSampleFilter< TSample >::DataObjectPointer
WeightedCovarianceSampleFilter< TSample >
::MakeOutput(DataObjectPointerArraySizeType index)
{
  MeasurementVectorSizeType measurementVectorSize = this->GetMeasurementVectorSize();

  if ( index == 0 )
    {
    MatrixType covarianceMatrix(measurementVectorSize, measurementVectorSize);
    covarianceMatrix.SetIdentity();
    typename MatrixDecoratedType::Pointer decoratedCovarianceMatrix = MatrixDecoratedType::New();
    decoratedCovarianceMatrix->Set(covarianceMatrix);
    return decoratedCovarianceMatrix.GetPointer();
    }

  if ( index >= 1 )
    {
    MeasurementVectorRealType mean;
    (void)mean; // for complainty pants : valgrind
    NumericTraits<MeasurementVectorRealType>::SetLength(mean, this->GetMeasurementVectorSize());
    // NumericTraits::SetLength also initializes array to zero
    typename MeasurementVectorDecoratedType::Pointer decoratedMean = MeasurementVectorDecoratedType::New();
    decoratedMean->Set( mean );
    return decoratedMean.GetPointer();
    }
  itkExceptionMacro("Trying to create output of index " << index << " larger than the number of output");
}

template< typename TSample >
inline void
WeightedCovarianceSampleFilter< TSample >
::ComputeCovarianceMatrixWithWeights()
{
  // set up input / output
  const SampleType *input = this->GetInput();

  MeasurementVectorSizeType measurementVectorSize = input->GetMeasurementVectorSize();

  MatrixDecoratedType *decoratedOutput =
    itkDynamicCastInDebugMode< MatrixDecoratedType * >( this->ProcessObject::GetOutput(0) );

  MatrixType output = decoratedOutput->Get();
  output.SetSize( measurementVectorSize, measurementVectorSize );
  output.Fill( NumericTraits< typename MatrixType::ValueType >::Zero );

  MeasurementVectorDecoratedType *decoratedMeanOutput =
    itkDynamicCastInDebugMode< MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(1) );

  WeightArrayType weightsArray(this->GetWeights());
  std::vector<std::vector<MeasurementType>> sampleComponents(measurementVectorSize);

  typename SampleType::ConstIterator iter =      input->Begin();
  const typename SampleType::ConstIterator end = input->End();

  MeasurementVectorType mean;
  MeasurementVectorType pbottom, ptop;
  itk::NumericTraits<MeasurementVectorType>::SetLength( pbottom, measurementVectorSize );
  itk::NumericTraits<MeasurementVectorType>::SetLength( ptop, measurementVectorSize );
  pbottom.Fill(itk::NumericTraits<MeasurementType>::min());
  ptop.Fill(itk::NumericTraits<MeasurementType>::max());

  if( this->m_RemoveOutliers ) {
	  MeasurementVectorType median;
	  // calculate percentiles
	  for ( unsigned int sampleVectorIndex = 0; iter != end; ++iter, ++sampleVectorIndex ) {
	    const MeasurementVectorType & measurement = iter.GetMeasurementVector();
	    const WeightValueType rawWeight = weightsArray[sampleVectorIndex];

	    if( rawWeight >= 0.9 ) {
	    	for(size_t c = 0; c < measurementVectorSize; c++)
	    		sampleComponents[c].push_back(measurement[c]);
	    }
	  }

	  itk::NumericTraits<MeasurementVectorType>::SetLength( median, measurementVectorSize );

	  size_t sampleSize = sampleComponents[0].size() - 1;
	  for(size_t c = 0; c < measurementVectorSize; c++) {
		  std::sort(sampleComponents[c].begin(), sampleComponents[c].end());
	  	  pbottom[c] = sampleComponents[c][int(0.02 * sampleSize)];
	  	  ptop[c] = sampleComponents[c][int(0.98 * sampleSize)];
	  	  median[c] = sampleComponents[c][int(0.50 * sampleSize)];
	  }
	  // mean = median;
  }

  MeasurementVectorDecoratedType *decoratedRangeMaxOutput =
	    itkDynamicCastInDebugMode< MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(2) );
  decoratedRangeMaxOutput->Set( ptop );

  MeasurementVectorDecoratedType *decoratedRangeMinOutput =
	    itkDynamicCastInDebugMode< MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(3) );
  decoratedRangeMinOutput->Set( pbottom );

  // reset sample iterator
  iter = input->Begin();
  // fills the lower triangle and the diagonal cells in the covariance matrix
  for ( unsigned int sampleVectorIndex = 0;
        iter != end;
        ++iter, ++sampleVectorIndex ) {
	  const MeasurementVectorType & measurement = iter.GetMeasurementVector();
	  WeightValueType rawWeight = weightsArray[sampleVectorIndex];

	  if (rawWeight > 0.0) {
		    for(size_t c = 0; c < measurementVectorSize; c++) {
		    	if (measurement[c] > ptop[c] || measurement[c] < pbottom[c] ) {
		    		weightsArray[sampleVectorIndex] = 0.0;
		    		continue;
		    	}
		    }
	  }
  }

  // calculate mean
  typedef itk::Statistics::WeightedMeanSampleFilter< SampleType > WeightedMeanFilterType;
  typename WeightedMeanFilterType::Pointer meanFilter = WeightedMeanFilterType::New();

  meanFilter->SetInput( input );
  meanFilter->SetWeights( weightsArray );
  meanFilter->Update();

  mean = meanFilter->GetMean();
  decoratedMeanOutput->Set( mean );

  // covariance algorithm
  MeasurementVectorRealType diff;
  itk::NumericTraits<MeasurementVectorRealType>::SetLength( diff, measurementVectorSize );

  WeightValueType totalWeight = itk::NumericTraits< WeightValueType >::Zero;
  WeightValueType totalSquaredWeight = itk::NumericTraits< WeightValueType >::Zero;

  // reset sample iterator
  iter = input->Begin();

  // fills the lower triangle and the diagonal cells in the covariance matrix
  for ( unsigned int sampleVectorIndex = 0;
        iter != end;
        ++iter, ++sampleVectorIndex )
    {
    const MeasurementVectorType & measurement = iter.GetMeasurementVector();
    const typename SampleType::AbsoluteFrequencyType frequency = iter.GetFrequency();
    const WeightValueType rawWeight = weightsArray[sampleVectorIndex];

    if ( rawWeight < 1.0e-8 )
    	continue;

    WeightValueType weight = ( rawWeight * static_cast< WeightValueType >( frequency ) );

    totalWeight += weight;
    totalSquaredWeight += ( weight * weight );

    for ( unsigned int dim = 0; dim < measurementVectorSize; ++dim )
      {
      const MeasurementRealType component =
        static_cast< MeasurementRealType >( measurement[dim] );

      diff[dim] = ( component - mean[dim] );
      }

    // updates the covariance matrix
    for ( unsigned int row = 0; row < measurementVectorSize; ++row )
      {
      for ( unsigned int col = 0; col < row + 1; ++col )
        {
        output(row, col) +=
          ( static_cast< MeasurementRealType >( weight ) * diff[row] * diff[col] );
        }
      }
    }

  // fills the upper triangle using the lower triangle
  for ( unsigned int row = 1; row < measurementVectorSize; ++row )
    {
    for ( unsigned int col = 0; col < row; ++col ) {
    	output(col, row) = output(row, col);
      }
    }

  const double normalizationFactor = ( totalWeight - ( totalSquaredWeight / totalWeight ) );

  if( normalizationFactor > vnl_math::eps )
    {
    const double inverseNormalizationFactor = 1.0 / normalizationFactor;

    output *= inverseNormalizationFactor;
    }
  else
    {
    itkExceptionMacro("Normalization factor was too close to zero. Value = " << normalizationFactor );
    }

  // Bias estimation
  // See http://en.wikipedia.org/wiki/Estimation_of_covariance_matrices#Bias_of_the_sample_covariance_matrix
  float p = measurementVectorSize;
  float n = totalWeight;
  float beta = (1/p) * (p * log(n) + p - bm::digamma(n-p+1) + (n - p + 1) * bm::digamma(n - p + 2) + bm::digamma(n+1) - (n+1)* bm::digamma(n+2));

  MatrixType ident;
  ident.SetSize( measurementVectorSize, measurementVectorSize );
  ident.SetIdentity();
  output+= output * exp( -beta );

  decoratedOutput->Set( output );
}

template< typename TSample >
const typename WeightedCovarianceSampleFilter< TSample >::MeasurementVectorDecoratedType *
WeightedCovarianceSampleFilter< TSample >
::GetRangeMaxOutput() const
{
  return static_cast< const MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(2) );
}

template< typename TSample >
const typename WeightedCovarianceSampleFilter< TSample >::MeasurementVectorRealType
WeightedCovarianceSampleFilter< TSample >
::GetRangeMax() const
{
  return this->GetRangeMaxOutput()->Get();
}

template< typename TSample >
const typename WeightedCovarianceSampleFilter< TSample >::MeasurementVectorDecoratedType *
WeightedCovarianceSampleFilter< TSample >
::GetRangeMinOutput() const
{
  return static_cast< const MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(3) );
}

template< typename TSample >
const typename WeightedCovarianceSampleFilter< TSample >::MeasurementVectorRealType
WeightedCovarianceSampleFilter< TSample >
::GetRangeMin() const
{
  return this->GetRangeMinOutput()->Get();
}

} // end of namespace itk

#endif
