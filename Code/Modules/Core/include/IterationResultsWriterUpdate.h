/*
 * IterationResultsWriterUpdate.h
 *
 *  Created on: Oct 16, 2013
 *      Author: oesteban
 */

#ifndef ITERATIONRESULTSWRITERUPDATE_H_
#define ITERATIONRESULTSWRITERUPDATE_H_

#include "IterationUpdate.h"

#include <itkImageFileWriter.h>
#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"



namespace rstk {

template< typename TOptimizer >
class IterationResultWriterUpdate: public IterationUpdate<TOptimizer>
{
public:
	typedef IterationResultWriterUpdate                 Self;
	typedef TOptimizer                                  OptimizerType;
	typedef typename itk::WeakPointer<OptimizerType>    OptimizerPointer;
	typedef IterationUpdate<TOptimizer>                 Superclass;
	typedef itk::SmartPointer<Self>                     Pointer;
	typedef itk::SmartPointer< const Self >             ConstPointer;

	typedef typename OptimizerType::ParametersType      ParametersType;
	typedef typename OptimizerType::FunctionalType      FunctionalType;
	typedef typename FunctionalType::ProbabilityMapType ProbabilityMapType;

	typedef rstk::DisplacementFieldComponentsFileWriter<ParametersType> ComponentsWriter;
	typedef itk::ImageFileWriter< ProbabilityMapType >  MapWriter;

	itkTypeMacro( IterationResultWriterUpdate, IterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(const itk::Object * object, const itk::EventObject & event) {
    	if( (typeid( event ) == typeid( itk::StartEvent ))  || (typeid( event ) == typeid( itk::IterationEvent ))) {
       		size_t nContours =this->m_Optimizer->GetFunctional()->GetCurrentContours().size();
        	std::stringstream ss;
       		for( size_t r = 0; r <= nContours; r++){
        		ss.str("");
        		ss << this->m_Prefix << "_region_" << r << "_it" << std::setfill('0')<<std::setw(3) << this->m_Optimizer->GetCurrentIteration() << ".nii.gz";
        		typename MapWriter::Pointer wr = MapWriter::New();
        		wr->SetInput( this->m_Optimizer->GetFunctional()->GetCurrentMap(r));
        		wr->SetFileName(ss.str().c_str() );
        		wr->Update();
        	}
    	}

    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {
        	typename ComponentsWriter::Pointer p = ComponentsWriter::New();
        	std::stringstream ss;
        	ss.str("");
        	ss << this->m_Prefix << "_parameters_" << std::setfill('0')  << std::setw(3) << this->m_Optimizer->GetCurrentIteration();
        	p->SetFileName( ss.str().c_str() );
        	p->SetInput( this->m_Optimizer->GetParameters() );
        	p->Update();

     		typename ComponentsWriter::Pointer f = ComponentsWriter::New();
    		ss.str("");
    		ss << this->m_Prefix << "_field_" << std::setfill('0')  << std::setw(3) << this->m_Optimizer->GetCurrentIteration();
    		f->SetFileName( ss.str().c_str() );
    		f->SetInput( this->m_Optimizer->GetCurrentDisplacementField() );
    		f->Update();
    	}
    }

    itkSetMacro( Prefix, std::string );
    itkGetConstMacro( Prefix, std::string );

    void SetOptimizer( OptimizerType * optimizer ) {
      m_Optimizer = optimizer;
      m_Optimizer->AddObserver( itk::IterationEvent(), this );
      m_Optimizer->AddObserver( itk::StartEvent(), this );
      m_Optimizer->AddObserver( itk::EndEvent(), this );
    }

protected:
    IterationResultWriterUpdate(){}

private:
    IterationResultWriterUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	OptimizerPointer   m_Optimizer;
	std::string        m_Prefix;
};

} // end namespace rstk



#endif /* ITERATIONRESULTSWRITERUPDATE_H_ */
