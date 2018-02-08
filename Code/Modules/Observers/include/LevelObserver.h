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

    void Execute(itk::Object *caller, const itk::EventObject & event) override {
        Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) override {
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
