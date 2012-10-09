/* --------------------------------------------------------------------------------------
 * File:    EnergyTerm.h
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

#ifndef EnergyTerm_h
#define EnergyTerm_h

#include <itkObject.h>

namespace rstk
{

template< class TDisplacementField, class TEnergyValue = float >
class EnergyTerm: public itk::Object {
public:
	typedef EnergyTerm                                     Self;
	typedef itk::Object                                    Superclass;
	typedef itk::SmartPointer<Self>                        Pointer;
	typedef itk::SmartPointer< const Self >                ConstPointer;
	typedef TEnergyValue                                   EnergyValueType;

	typedef TDisplacementField                             DisplacementFieldType;
	typedef typename DisplacementFieldType::Pointer        DisplacementFieldPointer;
	typedef typename DisplacementFieldType::ConstPointer   DisplacementFieldConstPointer;

	itkTypeMacro( EnergyTerm, itk::Object );




	itkSetConstObjectMacro( DisplacementField, DisplacementFieldType);
	itkGetObjectMacro( DisplacementField, DisplacementFieldType);

	/*
	 * Pure virtual method for energy term computation
	 */
	virtual EnergyValueType ComputeEnergy() = 0;

	virtual EnergyValueType GetGradient() = 0;


protected:
	EnergyTerm();
	virtual ~EnergyTerm() {}

	void PrintSelf( std::ostream &os, itk::Indent indent ) const {
		Superclass::PrintSelf( os, indent );

		// TODO PrintSelf
	}

	DisplacementFieldConstPointer m_DisplacementField;

private:
	EnergyTerm( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

}; // End of Class

} // End of namespace rstk

#endif /* EnergyTerm_h */
