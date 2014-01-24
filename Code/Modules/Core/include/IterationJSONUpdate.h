/*
 * IterationJSONUpdate.h
 *
 *  Created on: Oct 6, 2013
 *      Author: oesteban
 */

#ifndef ITERATIONJSONUPDATE_H_
#define ITERATIONJSONUPDATE_H_

#include "IterationUpdate.h"

#include <jsoncpp/json/json.h>

namespace rstk {

template< typename TOptimizer >
class IterationJSONUpdate: public IterationUpdate<TOptimizer>
{
public:
	typedef IterationJSONUpdate                       Self;
	typedef TOptimizer                                OptimizerType;
	typedef typename itk::WeakPointer<OptimizerType>  OptimizerPointer;
	typedef IterationUpdate<TOptimizer>               Superclass;
	typedef itk::SmartPointer<Self>                   Pointer;
	typedef itk::SmartPointer< const Self >           ConstPointer;
	typedef Json::Value                               JSONValue;

	itkTypeMacro( IterationJSONUpdate, IterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(const itk::Object * object, const itk::EventObject & event) {

		Json::Value iteration( Json::objectValue );

    	if( typeid( event ) == typeid( itk::StartEvent ) ) {
    		if( this->GetTrackEnergy() ) {
    			iteration["id"] = 0;
				iteration["energy"]["total"] = this->m_Optimizer->GetCurrentValue();
				iteration["energy"]["data"] = this->m_Optimizer->GetFunctional()->GetValue();
				iteration["energy"]["regularization"] = this->m_Optimizer->GetCurrentRegularizationEnergy();
				iteration["norm"] = 0.0;
				this->m_JSONRoot.append( iteration );
    		}
    	}

    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {

    		iteration["id"] = static_cast<Json::UInt64> (this->m_Optimizer->GetCurrentIteration());
    		if ( this->GetTrackEnergy() ) {
    			iteration["energy"]["total"] = this->m_Optimizer->GetCurrentEnergy();
    			iteration["energy"]["data"] = this->m_Optimizer->GetFunctional()->GetValue();
    			iteration["energy"]["regularization"] = this->m_Optimizer->GetCurrentRegularizationEnergy();
    		}
    		iteration["norm"] = this->m_Optimizer->GetCurrentValue();
    		this->m_JSONRoot.append( iteration );
    	}
    }

    // itkSetMacro( JSONRoot, JSONValue );
    itkGetConstMacro( JSONRoot, JSONValue );

    void SetOptimizer( OptimizerType * optimizer ) {
      m_Optimizer = optimizer;
      m_Optimizer->AddObserver( itk::IterationEvent(), this );
      m_Optimizer->AddObserver( itk::StartEvent(), this );
      m_Optimizer->AddObserver( itk::EndEvent(), this );
    }

protected:
    IterationJSONUpdate(){
    	m_JSONRoot = Json::Value( Json::arrayValue );
    }
    ~IterationJSONUpdate(){}

private:
    IterationJSONUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	JSONValue m_JSONRoot;
	OptimizerPointer   m_Optimizer;
};

} // end namespace rstk


#endif /* ITERATIONJSONUPDATE_H_ */
