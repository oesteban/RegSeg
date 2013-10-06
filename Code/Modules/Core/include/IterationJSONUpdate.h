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
    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {
    		Json::Value iteration( Json::objectValue );

    		iteration["id"] = static_cast<Json::UInt64> (this->m_Optimizer->GetCurrentIteration());
    		if ( this->GetTrackEnergy() ) {
    			iteration["energy"]["total"] = this->m_Optimizer->GetCurrentEnergy();
    			iteration["energy"]["functional"] = this->m_Optimizer->GetFunctional()->GetValue();
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
    }

protected:
    IterationJSONUpdate(){
    	m_JSONRoot = Json::Value( Json::arrayValue );
    }

private:
    IterationJSONUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	JSONValue m_JSONRoot;
	OptimizerPointer   m_Optimizer;
};

} // end namespace rstk


#endif /* ITERATIONJSONUPDATE_H_ */
