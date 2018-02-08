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

#ifndef SOURCE_DIRECTORY__MODULES_FUNCTIONAL_INCLUDE_UNIFORMMEMBERSHIPFUNCTION_H_
#define SOURCE_DIRECTORY__MODULES_FUNCTIONAL_INCLUDE_UNIFORMMEMBERSHIPFUNCTION_H_

#include <itkVariableSizeMatrix.h>
#include <itkMembershipFunctionBase.h>

namespace rstk {

template< typename TVector >
class UniformMembershipFunction:
  public itk::Statistics::MembershipFunctionBase< TVector >
{
public:
  /** Standard class typedefs */
  typedef UniformMembershipFunction Self;
  typedef itk::Statistics::MembershipFunctionBase< TVector >     Superclass;
  typedef itk::SmartPointer< Self >                              Pointer;
  typedef itk::SmartPointer< const Self >                        ConstPointer;

  /** Strandard macros */
  itkTypeMacro(UniformMembershipFunction, itk::Statistics::MembershipFunctionBase);
  itkNewMacro(Self);

  /** SmartPointer class for superclass */
  typedef typename Superclass::Pointer MembershipFunctionPointer;

  /** Typedef alias for the measurement vectors */
  typedef TVector MeasurementVectorType;

  /** Typedef to represent the length of measurement vectors */
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;

  /** Type of the mean vector. RealType on a vector-type is the same
   * vector-type but with a real element type.  */
  typedef typename itk::NumericTraits< MeasurementVectorType >::RealType MeasurementVectorRealType;
  typedef MeasurementVectorRealType  MeanVectorType;
  typedef typename itk::NumericTraits<MeasurementVectorRealType>::ScalarRealType MeasurementValueType;


  /**
   * Evaluate the Mahalanobis distance of a measurement using the
   * prescribed mean and covariance. Note that the Mahalanobis
   * distance is not a probability density. The square of the
   * distance is returned. */
  double Evaluate(const MeasurementVectorType & measurement) const override { return m_Value; }

  itkSetMacro(Value, double);
  itkGetMacro(Value, double);

  double GetOffsetTerm() const { return 0.0; }

protected:
  UniformMembershipFunction(void): m_Value(0.0), Superclass() {}
  virtual ~UniformMembershipFunction(void) {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const override { Superclass::PrintSelf(os, indent); }

private:
  double  m_Value;

};
} // end namespace itk


#endif /* SOURCE_DIRECTORY__MODULES_FUNCTIONAL_INCLUDE_UNIFORMMEMBERSHIPFUNCTION_H_ */
