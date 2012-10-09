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

#include <vector>
#include "EnergyTerm.h"
#include <itkObject.h>

namespace rstk
{

template< class TDisplacementField, class TEnergyValue = float >
class Energy: public EnergyTerm< TDisplacementField, TEnergyValue > {
public:
	typedef Energy                                    Self;
	typedef itk::rstk::EnergyTerm
			< TDisplacementField, TEnergyValue >      Superclass;
	typedef SmartPointer<Self>                        Pointer;
	typedef SmartPointer< const Self >                ConstPointer;

	itkTypeMacro( Energy, EnergyTerm );
	itkNewMacro( Self );

	typedef TDisplacementField                        DisplacementFieldType;
	typedef DisplacementFieldType::Pointer            DisplacementFieldPointer;
	typedef DisplacementFieldType::ConstPointer       DisplacementFieldConstPointer;

	typedef TEnergyValue                              EnergyValueType;

	typedef std::vector< const Self* >                EnergyComponentsContainer;

	EnergyValueType ComputeEnergy() {
		EnergyValueType result = 0.0;
		EnergyComponentsContainer::iterator end = this->m_EnergyComponentsContainer.end();

		for(EnergyComponentsContainer::iterator it = this->m_EnergyComponentsContainer.begin();
				it != end; it++ ) {
			result+= (*it)->ComputeEnergy();
		}
		return result;
	}

	EnergyValueType GetGradient() {
		float result = 0.0;
		EnergyComponentsContainer::iterator end = this->m_EnergyComponentsContainer.end();

		for(EnergyComponentsContainer::iterator it = this->m_EnergyComponentsContainer.begin();
				it != end; it++ ) {
			result+= (*it)->GetGradient();
		}
		return result;
	}


	void SetEnergyTerm( Superclass* energyTerm ) {
		energyTerm->SetDisplacementField( this->m_DisplacementField );
		this->m_EnergyComponentsContainer.push_back( static_cast< const Self* >( energyTerm ) );
	}


protected:
	Energy();
	virtual ~Energy() {}

	void PrintSelf( std::ostream &os, Indent indent ) const {
		Superclass::PrintSelf( os, indent );

		os << std::endl << "Energy Terms: " << std::flush;

		EnergyComponentsContainer::const_iterator end = this->m_EnergyComponentsContainer.end();
		for(EnergyComponentsContainer::const_iterator it = this->m_EnergyComponentsContainer.begin();
				it != end; it++ ) {
			os << std::endl << (*it)->PrintSelf(os, indent );
		}

		os << std::endl;
	}

private:
	EnergyComponentsContainer     m_EnergyComponentsContainer;

	Energy( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "Energy.txx"
#endif


#endif /* Energy_h */
