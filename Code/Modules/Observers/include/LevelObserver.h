// --------------------------------------------------------------------------------------
// File:          LevelObserver.h
// Date:          Jan 29, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
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
#include <itkMeshFileWriter.h>
#include "ACWERegistrationMethod.h"
#include "rstkVTKPolyDataWriter.h"


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

	typedef typename RegistrationMethodType::PriorsList        ContourList;
	typedef typename RegistrationMethodType::PriorsType        PriorsType;
	typedef rstk::VTKPolyDataWriter<PriorsType>                WriterType;
	// typedef itk::MeshFileWriter<PriorsType>                    WriterType;

	typedef typename RegistrationMethodType::TransformType       TransformType;
	typedef typename TransformType::AltCoeffType                 AltCoeffType;
	typedef rstk::CoefficientsWriter< AltCoeffType >             CoeffWriter;

	itkTypeMacro( LevelObserver, itk::CommandIterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(itk::Object *caller, const itk::EventObject & event) {
        Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) {
    	if (m_RegistrationMethod->GetVerbosity()==0) return;

    	std::stringstream ss;

    	if ( this->m_Prefix.size() > 0 ) {
    		char ch = *this->m_Prefix.rbegin();
    		if( ch != '_' )
    			this->m_Prefix.append( "_" );
    	}

    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {
    		// Contours and regions
    		ContourList conts = m_RegistrationMethod->GetCurrentContours();
    	    size_t nCont = conts.size();
    	    for ( size_t contid = 0; contid < nCont; contid++) {
    	    	typename WriterType::Pointer polyDataWriter = WriterType::New();
    	    	ss.str("");
    	    	ss << this->m_Prefix << "swarped_lev" << m_RegistrationMethod->GetCurrentLevel() << "_cont" << contid << ".vtk";
    	    	polyDataWriter->SetInput( conts[contid] );
    	    	polyDataWriter->SetFileName( ss.str().c_str() );
    	    	polyDataWriter->Update();
    	    }

    	    // Write transform parameters
    	    ss.str("");
    	    ss << this->m_Prefix << "coeff_" << m_RegistrationMethod->GetCurrentLevel() << ".vtu";
    	    typename CoeffWriter::Pointer w = CoeffWriter::New();
    	    w->SetFileName(ss.str().c_str());
    	    w->SetInput(m_RegistrationMethod->GetOptimizer()->GetTransform()->GetFlatParameters());
    	    w->Update();
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
