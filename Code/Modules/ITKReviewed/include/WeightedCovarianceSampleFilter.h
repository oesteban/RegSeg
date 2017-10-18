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

#ifndef __WeightedCovarianceSampleFilter_h
#define __WeightedCovarianceSampleFilter_h

#include <itkFunctionBase.h>
#include <itkCovarianceSampleFilter.h>
#include <itkDataObjectDecorator.h>

namespace itk
{
/** \class WeightedCovarianceSampleFilter
 * \brief Calculates the covariance matrix of the target sample data.
 *  where each measurement vector has an associated weight value
 *
 * Weight values can be specified in two ways: using a weighting function
 * or an array containing weight values. If none of these two is specified,
 * the covariance matrix is generated with equal weights.
 *
 * \sa CovarianceSampleFilter
 *
 * \ingroup ITKStatistics
 */

template< typename TSample >
class WeightedCovarianceSampleFilter:
		public itk::Statistics::CovarianceSampleFilter< TSample >
{
public:
  /** Standard class typedefs. */
  typedef WeightedCovarianceSampleFilter         Self;
  typedef itk::Statistics::CovarianceSampleFilter< TSample >
                                                 Superclass;
  typedef itk::SmartPointer< Self >              Pointer;
  typedef itk::SmartPointer< const Self >        ConstPointer;

  /** Standard Macros */
  itkTypeMacro(WeightedCovarianceSampleFilter, itk::Statistics::CovarianceSampleFilter);
  itkNewMacro(Self);

  /** Types derived from the base class */
  typedef typename Superclass::SampleType                     SampleType;
  typedef typename Superclass::MeasurementVectorType          MeasurementVectorType;
  typedef typename Superclass::MeasurementVectorSizeType      MeasurementVectorSizeType;
  typedef typename Superclass::MeasurementType                MeasurementType;

  /** Types derived from the base class */
  typedef typename Superclass::MeasurementVectorRealType      MeasurementVectorRealType;
  typedef typename Superclass::MeasurementRealType            MeasurementRealType;


  /** Type of weight values */
  typedef double WeightValueType;


  /** Array type for weights */
  typedef itk::Array< WeightValueType > WeightArrayType;

  /** Type of DataObjects to use for the weight array type */
  typedef itk::SimpleDataObjectDecorator< WeightArrayType > InputWeightArrayObjectType;

  /** Method to set the input value of the weight array */
  itkSetGetDecoratedInputMacro(Weights, WeightArrayType);


  /** Weight calculation function type */
  typedef itk::FunctionBase< MeasurementVectorType, WeightValueType > WeightingFunctionType;

  /** Type of DataObjects to use for Weight function */
  typedef itk::DataObjectDecorator< WeightingFunctionType > InputWeightingFunctionObjectType;

  /** Method to set/get the weighting function */
  itkSetGetDecoratedObjectInputMacro(WeightingFunction, WeightingFunctionType);

  /** Types derived from the base class */
  typedef typename Superclass::MatrixType          MatrixType;
  typedef typename Superclass::MatrixDecoratedType MatrixDecoratedType;

  /** Types derived from the base class */
  typedef typename Superclass::MeasurementVectorDecoratedType MeasurementVectorDecoratedType;
  typedef typename Superclass::OutputType                     OutputType;

  const MeasurementVectorRealType GetRangeMax() const;
  const MeasurementVectorDecoratedType * GetRangeMaxOutput() const;
  const MeasurementVectorRealType GetRangeMin() const;
  const MeasurementVectorDecoratedType * GetRangeMinOutput() const;

protected:
  WeightedCovarianceSampleFilter();
  virtual ~WeightedCovarianceSampleFilter();
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

  /** DataObject pointer */
  typedef itk::DataObject::Pointer DataObjectPointer;

  typedef itk::ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

  void GenerateData();

  /** Compute covariance matrix with weights computed from a function */
  void ComputeCovarianceMatrixWithWeightingFunction();

  /** Compute covariance matrix with weights specified in an array */
  void ComputeCovarianceMatrixWithWeights();

private:
  WeightedCovarianceSampleFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                 //purposely not implemented

  bool m_RemoveOutliers;

};  // end of class
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "WeightedCovarianceSampleFilter.hxx"
#endif

#endif
