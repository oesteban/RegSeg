// --------------------------------------------------------------------------------------
// File:          LevelObserver.h
// Date:          Jan 29, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2015, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWEReg
//
// ACWEReg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWEReg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWEReg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef SOURCE_DIRECTORY__MODULES_OBSERVERS_INCLUDE_LEVELOBSERVER_H_
#define SOURCE_DIRECTORY__MODULES_OBSERVERS_INCLUDE_LEVELOBSERVER_H_

#include <itkCommand.h>
#include <itkWeakPointer.h>
#include "ACWERegistrationMethod.h"


namespace rstk {

template< typename TRegistrationMethod >
class LevelObserver: public itk::Command
{
public:
	typedef LevelObserver                                      Self;
	typedef TRegistrationMethod                                RegistrationMethodType;
	typedef typename itk::WeakPointer<RegistrationMethodType>  RegistrationMethodPointer;
	typedef itk::Command                                       Superclass;
	typedef itk::SmartPointer<Self>                            Pointer;
	typedef itk::SmartPointer< const Self >                    ConstPointer;

	itkTypeMacro( LevelObserver, itk::CommandIterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(itk::Object *caller, const itk::EventObject & event) {
        Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) {
    	std::stringstream ss;

    	if ( this->m_Prefix.size() > 0 ) {
    		char ch = *this->m_Prefix.rbegin();
    		if( ch != '_' )
    			this->m_Prefix.append( "_" );
    	}

    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {
    		ss.str("");
    		ss << this->m_Prefix << "uk_" << std::setfill('0') << "lev" << this->m_RegistrationMethod->GetCurrentLevel();

    		std::cout << "Here I'd write ->" << ss.str() << std::endl;

    	}
    }

    void SetRegistrationMethod( RegistrationMethodType * _arg ) {
      m_RegistrationMethod = _arg;
      m_RegistrationMethod->AddObserver( itk::IterationEvent(), this );
    }

    itkSetMacro( Prefix, std::string );
    itkGetConstMacro( Prefix, std::string );
protected:
	LevelObserver(): m_Prefix("") {}
	~LevelObserver(){}

private:
	LevelObserver( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented
	RegistrationMethodPointer   m_RegistrationMethod;
	std::string                 m_Prefix;
};

} // end namespace rstk


#endif /* SOURCE_DIRECTORY__MODULES_OBSERVERS_INCLUDE_LEVELOBSERVER_H_ */
