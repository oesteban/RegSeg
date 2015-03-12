// --------------------------------------------------------------------------------------
// File:          BSplineSecondDerivativeKernelFunction.h
// Date:          Feb 6, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

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
