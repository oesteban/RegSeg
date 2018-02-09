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

#ifndef ITERATIONSTDOUTUPDATE_H_
#define ITERATIONSTDOUTUPDATE_H_

#include "IterationUpdate.h"

#include <boost/lexical_cast.hpp>
#include <ctime>

namespace rstk {

template< typename TOptimizer >
class IterationStdOutUpdate: public IterationUpdate<TOptimizer>
{
public:
	typedef IterationStdOutUpdate                       Self;
	typedef TOptimizer                                OptimizerType;
	typedef typename itk::WeakPointer<OptimizerType>  OptimizerPointer;
	typedef IterationUpdate<TOptimizer>               Superclass;
	typedef itk::SmartPointer<Self>                   Pointer;
	typedef itk::SmartPointer< const Self >           ConstPointer;

	itkTypeMacro( IterationStdOutUpdate, IterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(const itk::Object * object, const itk::EventObject & event) override {
		size_t it = this->m_Optimizer->GetCurrentIteration();

    	if( typeid( event ) == typeid( itk::StartEvent ) ) {
    		m_StartTime = clock();

    		std::cout << "[ini] " << std::setw(8) << "     N/A";
    		if( !this->m_Optimizer->GetUseLightWeightConvergenceChecking() ) {
    			std::cout << " " << this->m_Optimizer->GetCurrentEnergy();
				std::cout << " " << this->m_Optimizer->GetFunctional()->GetValue();
				std::cout << " " << this->m_Optimizer->GetCurrentRegularizationEnergy();
				std::cout << " ||";

				typename OptimizerType::FunctionalType::MeasureArray es = this->m_Optimizer->GetFunctional()->GetRegionValue();
				for (size_t r = 0; r < es.Size(); r++) {
					std::cout << " " << es[r];
				}
    		}
    		std::cout << "." << std::endl;
    	}

		if( typeid( event ) == typeid( itk::IterationEvent ) ) {
			std::cout << "[" << std::setw(3) << it << "] " << std::setw(8) << this->m_Optimizer->GetCurrentValue();
			std::cout << " " << this->m_Optimizer->GetMaximumGradient();
			std::cout << " " << this->m_Optimizer->GetConvergenceValue();
			std::cout << " ||";
			if( !this->m_Optimizer->GetUseLightWeightConvergenceChecking() ) {
    			std::cout << " " << this->m_Optimizer->GetCurrentEnergy();
				std::cout << " " << this->m_Optimizer->GetFunctional()->GetValue();
				std::cout << " " << this->m_Optimizer->GetCurrentRegularizationEnergy();
				std::cout << " ||";

				typename OptimizerType::FunctionalType::MeasureArray es = this->m_Optimizer->GetFunctional()->GetRegionValue();
				for (size_t r = 0; r < es.Size(); r++) {
					std::cout << " " << es[r];
				}

				std::cout << " ||";
    		}

			const std::vector< size_t > off = this->m_Optimizer->GetFunctional()->GetOffMaskVertices();
			std::cout << " OffMask=(";
			for (size_t c = 0; c<off.size(); c++) {
				std::cout << off[c] << ", ";
			}
			std::cout << ")";

    		std::cout << "." << std::endl;
		}

		if( typeid( event ) == typeid( FunctionalModifiedEvent ) )  {
			std::cout << "[" << std::setw(3) << it << "] Descriptors updated." << std::endl;
		}

		if( typeid( event ) == typeid( itk::EndEvent ) ) {
			std::cout << "[end] " << std::setw(8) << this->m_Optimizer->GetCurrentValue();
			if( !this->m_Optimizer->GetUseLightWeightConvergenceChecking() ) {
				std::cout << " " << this->m_Optimizer->GetCurrentEnergy();
				std::cout << " " << this->m_Optimizer->GetFunctional()->GetValue();
				std::cout << " " << this->m_Optimizer->GetCurrentRegularizationEnergy();
				std::cout << " ||";

				typename OptimizerType::FunctionalType::MeasureArray es = this->m_Optimizer->GetFunctional()->GetRegionValue();
				for (size_t r = 0; r < es.Size(); r++) {
					std::cout << " " << es[r];
				}
			}
			std::cout << "." << std::endl;

			m_StopTime = clock();
			float tot = (float) (((double) (m_StopTime - m_StartTime)) / CLOCKS_PER_SEC);
			std::cout << "Elapsed time: " << tot << "s." << std::endl;
		}
    }

    void SetOptimizer( OptimizerType * optimizer ) {
      m_Optimizer = optimizer;
      m_Optimizer->AddObserver( itk::IterationEvent(), this );
      m_Optimizer->AddObserver( itk::StartEvent(), this );
      m_Optimizer->AddObserver( itk::EndEvent(), this );
      m_Optimizer->AddObserver( itk::ModifiedEvent(), this );
      m_Optimizer->AddObserver( FunctionalModifiedEvent(), this );
    }

protected:
    IterationStdOutUpdate(): m_StartTime(0), m_StopTime(0) { }
    ~IterationStdOutUpdate(){}

private:
    IterationStdOutUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	OptimizerPointer   m_Optimizer;
	clock_t m_StartTime;
	clock_t m_StopTime;
};

} // end namespace rstk



#endif /* ITERATIONSTDOUTUPDATE_H_ */
