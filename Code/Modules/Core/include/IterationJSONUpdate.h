/*
 * IterationJSONUpdate.h
 *
 *  Created on: Oct 6, 2013
 *      Author: oesteban
 */

#ifndef ITERATIONJSONUPDATE_H_
#define ITERATIONJSONUPDATE_H_

#include "IterationUpdate.h"

#include <boost/lexical_cast.hpp>
#include <jsoncpp/json/json.h>
#include <ctime>

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
	typedef typename OptimizerType::InternalComputationValueType InternalOptimizerValue;

	itkTypeMacro( IterationJSONUpdate, IterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(const itk::Object * object, const itk::EventObject & event) {
		size_t it = this->m_Optimizer->GetCurrentIteration();
    	Json::Value itnode;

    	if( it > 0 && it==m_LastIt ) {
    		itnode = m_Last;
    	} else {
    		if ( it > 0 ) {
    			this->m_JSONRoot.append( m_Last );
    		}
    		itnode = Json::Value( Json::objectValue );
    		itnode["id"] = static_cast<Json::UInt64>( it );
    	}

    	if( typeid( event ) == typeid( itk::StartEvent ) ) {
    		m_StartTime = clock();

    		if( !this->m_Optimizer->GetUseLightWeightConvergenceChecking() ) {
				itnode["energy"]["total"] = this->m_Optimizer->GetCurrentEnergy();
				itnode["energy"]["data"] = this->m_Optimizer->GetFunctional()->GetValue();
				itnode["energy"]["regularization"] = this->m_Optimizer->GetCurrentRegularizationEnergy();

				typename OptimizerType::FunctionalType::MeasureArray es = this->m_Optimizer->GetFunctional()->GetRegionValue();
				JSONValue enode = Json::Value( Json::arrayValue );
				for (size_t r = 0; r < es.Size(); r++) {
					enode.append(Json::Value(es[r]));
				}
				itnode["energy"]["region"] = enode;
    		}
    		itnode["descriptors"] = this->ParseTree( this->m_Optimizer->GetFunctional()->PrintFormattedDescriptors() );
    		itnode["step_size"] = this->m_Optimizer->GetStepSize();

    		JSONValue size = Json::Value( Json::arrayValue );
    		JSONValue spacing = Json::Value( Json::arrayValue );
    		JSONValue maxdisp = Json::Value( Json::arrayValue );

    		for( size_t i = 0; i < this->m_Optimizer->GetGridSize().GetSizeDimension(); i++) {
    			size.append(  static_cast<Json::UInt64>( this->m_Optimizer->GetGridSize()[i] ) );
    			spacing.append(  static_cast<Json::UInt64>( this->m_Optimizer->GetGridSpacing()[i] ) );
    			maxdisp.append(  static_cast<Json::UInt64>( this->m_Optimizer->GetMaxDisplacement()[i] ) );
    		}

    		itnode["transform"]["grid-size"] = size;
    		itnode["transform"]["grid-spacing"] = spacing;
    		itnode["transform"]["max-displacement"] = maxdisp;
    		itnode["transform"]["num-threads"] = this->m_Optimizer->GetTransform()->GetNumberOfThreads();
    	}

		if( typeid( event ) == typeid( itk::IterationEvent ) ) {
			if( !this->m_Optimizer->GetUseLightWeightConvergenceChecking() ) {
				itnode["energy"]["total"] = this->m_Optimizer->GetCurrentEnergy();
				itnode["energy"]["data"] = this->m_Optimizer->GetFunctional()->GetValue();
				itnode["energy"]["regularization"] = this->m_Optimizer->GetCurrentRegularizationEnergy();

				typename OptimizerType::FunctionalType::MeasureArray es = this->m_Optimizer->GetFunctional()->GetRegionValue();
				JSONValue enode = Json::Value( Json::arrayValue );
				for (size_t r = 0; r < es.Size(); r++) {
					enode.append(Json::Value(es[r]));
				}
				itnode["energy"]["region"] = enode;
			}

			InternalOptimizerValue val = this->m_Optimizer->GetConvergenceValue();
			if( val < itk::NumericTraits<InternalOptimizerValue>::max() ) {
				itnode["convergence"]["value"] = val;
			} else {
				itnode["convergence"]["value"] = "inf";
			}
			itnode["convergence"]["norm"] = this->m_Optimizer->GetCurrentNorm();
			itnode["convergence"]["step_size"] = this->m_Optimizer->GetStepSize();
			itnode["convergence"]["max_gradient"] = this->m_Optimizer->GetMaximumGradient();
			itnode["convergence"]["momentum"] = this->m_Optimizer->GetMomentum();
			itnode["speed"]["max"] = this->m_Optimizer->GetMaxSpeed();
			itnode["speed"]["median"] = this->m_Optimizer->GetMeanSpeed();
			itnode["speed"]["average"] = this->m_Optimizer->GetAvgSpeed();


			std::vector< size_t > off = this->m_Optimizer->GetFunctional()->GetOffMaskVertices();
			JSONValue offnode = Json::Value( Json::arrayValue );
			for (size_t c = 0; c<off.size(); c++) {
				offnode.append(Json::UInt64(off[c]));
			}
			itnode["off-grid"] = offnode;

			itnode["diffemorphic"] = this->m_Optimizer->GetIsDiffeomorphic();
			itnode["diffemorphism-forced"] = this->m_Optimizer->GetDiffeomorphismForced();
		}

		if( typeid( event ) == typeid( FunctionalModifiedEvent ) )  {
			itnode["descriptors"] = this->ParseTree( this->m_Optimizer->GetFunctional()->PrintFormattedDescriptors() );
		}

		if( typeid( event ) == typeid( itk::EndEvent ) ) {
			itnode["convergence"]["norm"] = this->m_Optimizer->GetCurrentNorm();
			itnode["convergence"]["step_size"] = this->m_Optimizer->GetStepSize();
			itnode["convergence"]["max_gradient"] = this->m_Optimizer->GetMaximumGradient();
			itnode["convergence"]["momentum"] = this->m_Optimizer->GetMomentum();
			itnode["speed"]["max"] = this->m_Optimizer->GetMaxSpeed();
			itnode["speed"]["median"] = this->m_Optimizer->GetMeanSpeed();
			itnode["speed"]["average"] = this->m_Optimizer->GetAvgSpeed();
			m_StopTime = clock();
			float tot_t = (float) (((double) (m_StopTime - m_StartTime)) / CLOCKS_PER_SEC);
			// JSON Summary
			itnode["time"]["processing"] = tot_t;
			itnode["summary"]["energy"]["total"] = this->m_Optimizer->GetCurrentEnergy();
			itnode["summary"]["energy"]["data"] = this->m_Optimizer->GetFunctional()->GetValue();
			itnode["summary"]["energy"]["regularization"] = this->m_Optimizer->GetCurrentRegularizationEnergy();

			typename OptimizerType::FunctionalType::MeasureArray es = this->m_Optimizer->GetFunctional()->GetRegionValue();
			JSONValue enode = Json::Value( Json::arrayValue );
			for (size_t r = 0; r < es.Size(); r++) {
				enode.append(Json::Value(es[r]));
			}
			itnode["summary"]["energy"]["region"] = enode;

			itnode["summary"]["iterations"] = Json::Int (this->m_Optimizer->GetCurrentIteration());
			itnode["summary"]["conv_status"] = this->m_Optimizer->GetStopCondition();
			itnode["summary"]["stop_msg"] = this->m_Optimizer->GetStopConditionDescription();
			itnode["summary"]["is-diffeomorphic"] = Json::Int( this->m_Optimizer->GetIsDiffeomorphic() );
			this->m_JSONRoot.append( itnode );
		}

    	m_Last = itnode;
    	m_LastIt = it;
    }

    // itkSetMacro( JSONRoot, JSONValue );
    itkGetConstMacro( JSONRoot, JSONValue );

    void SetOptimizer( OptimizerType * optimizer ) {
      m_Optimizer = optimizer;
      m_Optimizer->AddObserver( itk::IterationEvent(), this );
      m_Optimizer->AddObserver( itk::StartEvent(), this );
      m_Optimizer->AddObserver( itk::EndEvent(), this );
      m_Optimizer->AddObserver( itk::ModifiedEvent(), this );
      m_Optimizer->AddObserver( FunctionalModifiedEvent(), this );
    }

protected:
    IterationJSONUpdate(): m_LastIt(0), m_StartTime(0), m_StopTime(0) {
    	m_JSONRoot = Json::Value( Json::arrayValue );
    	m_Last = Json::Value( Json::objectValue );
    }
    ~IterationJSONUpdate(){}

private:
    IterationJSONUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	Json::Value ParseTree( std::string str ){
		Json::Value node( Json::objectValue );
		Json::Value val( Json::arrayValue );
		Json::Reader r;

		if( r.parse( str, node, false ) ) {
			if ( node.size() > 0 ) {
				for( Json::ValueIterator itr = node.begin() ; itr != node.end() ; itr++ ) {
					val.append( *itr );
				}
			}
		}

		return val;
	}

	JSONValue m_JSONRoot;
	JSONValue m_Last;
	OptimizerPointer   m_Optimizer;
	size_t m_LastIt;
	clock_t m_StartTime;
	clock_t m_StopTime;
};

} // end namespace rstk


#endif /* ITERATIONJSONUPDATE_H_ */
