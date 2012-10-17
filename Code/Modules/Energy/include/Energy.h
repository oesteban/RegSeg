/* --------------------------------------------------------------------------------------
 * File:    Energy.h
 * Date:    19/02/2012
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

#ifndef Energy_h
#define Energy_h


#include "EnergyBase.h"
#include <vector>
#include <itkObject.h>

namespace rstk
{

class Energy: public EnergyBase {
public:
	typedef Energy                                    Self;
	typedef EnergyBase                                Superclass;
	typedef itk::SmartPointer<Self>                   Pointer;
	typedef itk::SmartPointer< const Self >           ConstPointer;

	itkTypeMacro( Energy, EnergyBase );
	itkNewMacro( Self );

	typedef Superclass::MeasureType                   MeasureType;
	typedef std::vector< const Superclass* >          EnergyTermsContainer;

	void SetEnergyTerm( Superclass* term);

	/** Initialize the Energy by making sure that all the components
	 *  are present and plugged together correctly, and initializing
	 *  internal variables as required. This is for one-time initialization,
	 *  e.g. before starting an optimization process. */
	virtual void Initialize(void) throw ( itk::ExceptionObject ) = 0;

	/** Type to represent the number of parameters that are being optimized at
	 * any given iteration of the optimizer. */
	typedef unsigned int NumberOfParametersType;

	/** Calculate and return the value for the metric based on the current
	 * transformation(s). The result is both returned, and stored in the
	 * m_Value member variable. */
	MeasureType GetValue() const;

	/**
	 * This method returns the derivative based on the current
	 * transformation(s). */
	virtual void GetDerivative( DerivativeType & ) const = 0;

	/** This method returns the derivative and value based on the current
	 * transformation(s). */
	virtual void GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const = 0;

	/** Methods for working with the metric's 'active' transform, e.g. the
	 * transform being optimized in the case of registration. Some of these are
	 * used in non-metric classes, e.g. optimizers. */
	virtual NumberOfParametersType GetNumberOfParameters() const = 0;
	virtual NumberOfParametersType GetNumberOfLocalParameters() const = 0;

	/** Set the active transform's parameters */
	virtual void SetParameters( ParametersType & params ) = 0;

	/** Get a const reference to the active transform's parameters */
	virtual const ParametersType & GetParameters() const = 0;

	/** Return whether the metric's active transform has local support,
	 * e.g. whether it is dense/high-dimensional. */
	virtual bool HasLocalSupport() const = 0;

	/** Update the parameters of the metric's active transform.
	 * Typically this call is passed through directly to the transform.
	 * \c factor is a scalar multiplier for each value in update, and
	 * defaults to 1.0 .
	 * \c derivative must be the proper size, as retrieved
	 * from GetNumberOfParameters. */
	virtual void UpdateTransformParameters( const DerivativeType & derivative,
			ParametersValueType factor = itk::NumericTraits<ParametersValueType>::One) = 0;


protected:
	Energy();
	virtual ~Energy() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const;

	ParametersType m_Parameters;

private:
	EnergyTermsContainer     m_EnergyTermsContainer;

	Energy( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk


#endif /* Energy_h */
