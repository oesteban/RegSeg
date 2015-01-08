/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __WeightedCovarianceSampleFilter_hxx
#define __WeightedCovarianceSampleFilter_hxx

#include "WeightedCovarianceSampleFilter.h"

namespace rstk
{
template< typename TSample >
WeightedCovarianceSampleFilter< TSample >
::WeightedCovarianceSampleFilter():
 m_MeanSet(false),
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

  const WeightArrayType & weightsArray = this->GetWeights();
  std::vector<MeasurementType> sampleComponents[measurementVectorSize];

  typename SampleType::ConstIterator iter =      input->Begin();
  const typename SampleType::ConstIterator end = input->End();

  MeasurementVectorType pbottom, ptop;
  itk::NumericTraits<MeasurementVectorType>::SetLength( pbottom, measurementVectorSize );
  itk::NumericTraits<MeasurementVectorType>::SetLength( ptop, measurementVectorSize );
  pbottom.Fill(itk::NumericTraits<MeasurementType>::min());
  ptop.Fill(itk::NumericTraits<MeasurementType>::max());


  MeasurementVectorType mean;
  if( this->m_MeanSet ) {
	  mean = decoratedMeanOutput->Get();
  } else {
	  // calculate mean
	  for ( unsigned int sampleVectorIndex = 0; iter != end; ++iter, ++sampleVectorIndex ) {
	    const MeasurementVectorType & measurement = iter.GetMeasurementVector();
	    const WeightValueType rawWeight = weightsArray[sampleVectorIndex];

	    if( rawWeight >= 0.9 ) {
	    	for(size_t c = 0; c < measurementVectorSize; c++)
	    		sampleComponents[c].push_back(measurement[c]);
	    }
	  }

	  itk::NumericTraits<MeasurementVectorType>::SetLength( mean, measurementVectorSize );

	  size_t sampleSize = sampleComponents[0].size() - 1;
	  for(size_t c = 0; c < measurementVectorSize; c++) {
		  std::sort(sampleComponents[c].begin(), sampleComponents[c].end());
	  	  pbottom[c] = sampleComponents[c][int(0.02 * sampleSize)];
	  	  ptop[c] = sampleComponents[c][int(0.98 * sampleSize)];
	  	  mean[c] = sampleComponents[c][int(0.50 * sampleSize)];
	  }
	  decoratedMeanOutput->Set( mean );
  }

  MeasurementVectorDecoratedType *decoratedRangeMaxOutput =
	    itkDynamicCastInDebugMode< MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(2) );
  decoratedRangeMaxOutput->Set( ptop );

  MeasurementVectorDecoratedType *decoratedRangeMinOutput =
	    itkDynamicCastInDebugMode< MeasurementVectorDecoratedType * >( this->ProcessObject::GetOutput(3) );
  decoratedRangeMinOutput->Set( pbottom );

  // covariance algorithm
  MeasurementVectorRealType diff;
  itk::NumericTraits<MeasurementVectorRealType>::SetLength( diff, measurementVectorSize );

  WeightValueType totalWeight = itk::NumericTraits< WeightValueType >::Zero;
  WeightValueType totalSquaredWeight = itk::NumericTraits< WeightValueType >::Zero;

  // reset sample iterator
  iter = input->Begin();

  bool isOutlier;
  // fills the lower triangle and the diagonal cells in the covariance matrix
  for ( unsigned int sampleVectorIndex = 0;
        iter != end;
        ++iter, ++sampleVectorIndex )
    {
    const MeasurementVectorType & measurement = iter.GetMeasurementVector();
    isOutlier = false;

    const typename SampleType::AbsoluteFrequencyType frequency = iter.GetFrequency();

    if (frequency != 1.0) {
    	std::cout << frequency << std::endl;
    }


    const WeightValueType rawWeight = weightsArray[sampleVectorIndex];

    WeightValueType weight = ( rawWeight * static_cast< WeightValueType >( frequency ) );

    for(size_t c = 0; c < measurementVectorSize; c++) {
    	if (measurement[c] > ptop[c] || measurement[c] < pbottom[c] ) {
    		weight = 0.0;
    		continue;
    	}
    }

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
