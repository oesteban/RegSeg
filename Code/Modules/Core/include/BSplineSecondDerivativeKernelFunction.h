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

#ifndef BSPLINESECONDDERIVATIVEKERNELFUNCTION_H_
#define BSPLINESECONDDERIVATIVEKERNELFUNCTION_H_

#include "itkBSplineKernelFunction.h"

namespace itk
{
/** \class BSplineSecondDerivativeKernelFunction
 * \brief Second Derivative of a BSpline kernel used for density estimation and
 *  nonparameteric regression.
 *
 * This class encapsulates the derivative of a BSpline kernel for
 * density estimation or nonparameteric regression.
 * See documentation for KernelFunctionBase for more details.
 *
 * This class is templated over the spline order.
 * \warning Evaluate is only implemented for spline order 1 to 4
 *
 * \sa KernelFunctionBase
 *
 * \ingroup Functions
 * \ingroup ITKCommon
 */
template< unsigned int VSplineOrder = 3, typename TRealValueType = double >
class ITK_EXPORT BSplineSecondDerivativeKernelFunction:public KernelFunctionBase<TRealValueType>
{
public:
  /** Standard class typedefs. */
  typedef BSplineSecondDerivativeKernelFunction    Self;
  typedef KernelFunctionBase<TRealValueType> Superclass;
  typedef SmartPointer< Self >               Pointer;

  typedef typename Superclass::RealType  RealType;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BSplineSecondDerivativeKernelFunction, KernelFunctionBase);

  /** Enum of for spline order. */
  itkStaticConstMacro(SplineOrder, unsigned int, VSplineOrder);

  /** Evaluate the function. */
  inline TRealValueType Evaluate( const TRealValueType & u ) const
    {
    return this->Evaluate( Dispatch< VSplineOrder >(), u );
    }

protected:
  BSplineSecondDerivativeKernelFunction() {}
  virtual ~BSplineSecondDerivativeKernelFunction(){}

  void PrintSelf(std::ostream & os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    os << indent  << "Spline Order: " << SplineOrder << std::endl;
    }

private:
  BSplineSecondDerivativeKernelFunction(const Self &); //purposely not implemented
  void operator=(const Self &);                  //purposely not implemented

  /** Structures to control overloaded versions of Evaluate */
  struct DispatchBase {};
  template< unsigned int >
  struct Dispatch: public DispatchBase {};

  /** Evaluate the function:  zeroth order spline. */
  inline TRealValueType Evaluate( const Dispatch<0>&, const TRealValueType & itkNotUsed( u ) )
    const
    {
    return NumericTraits< TRealValueType >::Zero;
    }

  /** Evaluate the function:  first order spline */
  inline TRealValueType Evaluate( const Dispatch<1>&, const TRealValueType & itkNotUsed( u ) ) const
    {
	return NumericTraits< TRealValueType >::Zero;
    }

  /** Evaluate the function:  second order spline. */
  inline TRealValueType Evaluate( const Dispatch<2>&, const TRealValueType & itkNotUsed( u ) ) const
    {
	  return NumericTraits< TRealValueType >::Zero;
    }

  /** Evaluate the function:  third order spline. */
  inline TRealValueType Evaluate( const Dispatch<3>&, const TRealValueType& u ) const
    {
    if( ( u > -NumericTraits< TRealValueType >::One ) && ( u < NumericTraits< TRealValueType >::One ) )
      {
      return ( static_cast< TRealValueType >(-2.0) + static_cast< TRealValueType >(3.0) * u);
      }
    else if( ( u >= NumericTraits< TRealValueType >::One ) && ( u < static_cast< TRealValueType >(2.0) ) )
      {
      return ( static_cast< TRealValueType >(2.0) + static_cast< TRealValueType >(-1.0) * u );
      }
    else if( ( u > static_cast< TRealValueType >(-2.0) ) && ( u <= -NumericTraits< TRealValueType >::One ) )
      {
      return ( static_cast< TRealValueType >(2.0) + static_cast< TRealValueType >(-1.0) * u );
      }
    else
      {
      return NumericTraits< TRealValueType >::Zero;
      }
    }

  /** Evaluate the function:  unimplemented spline order */
  inline TRealValueType Evaluate( const DispatchBase&, const TRealValueType& ) const
    {
    itkExceptionMacro( "Evaluate not implemented for spline order " << SplineOrder );
    return NumericTraits< TRealValueType >::Zero; // This is to avoid compiler warning about missing
    // return statement. It should never be evaluated.
    }
};
} // end namespace itk


#endif /* BSPLINESECONDDERIVATIVEKERNELFUNCTION_H_ */
