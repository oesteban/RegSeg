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

#ifndef ITERATIONUPDATE_H_
#define ITERATIONUPDATE_H_

#include <itkCommandIterationUpdate.h>
#include <itkWeakPointer.h>
#include "OptimizerBase.h"


namespace rstk {

template< typename TOptimizer >
class IterationUpdate: public itk::CommandIterationUpdate<TOptimizer>
{
public:
	typedef IterationUpdate                           Self;
	typedef TOptimizer                                OptimizerType;
	typedef typename itk::WeakPointer<OptimizerType>  OptimizerPointer;
	typedef itk::CommandIterationUpdate<TOptimizer>   Superclass;
	typedef itk::SmartPointer<Self>                   Pointer;
	typedef itk::SmartPointer< const Self >           ConstPointer;

	itkTypeMacro( IterationUpdate, itk::CommandIterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

    void Execute(itk::Object *caller, const itk::EventObject & event) override {
        Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) override {

    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {
    		std::cout << "[" << this->m_Optimizer->GetCurrentIteration() << "] ";
    		if( !this->m_Optimizer->GetUseLightWeightConvergenceChecking() ) {
    			std::cout << this->m_Optimizer->GetCurrentEnergy() << " | " << this->m_Optimizer->GetFunctional()->GetValue() << " | " << this->m_Optimizer->GetCurrentRegularizationEnergy();
    		}
    		std::cout << "\t" << this->m_Optimizer->GetCurrentValue() << std::endl;
    	}
    }

    void SetOptimizer( OptimizerType * optimizer ) {
      m_Optimizer = optimizer;
      m_Optimizer->AddObserver( itk::IterationEvent(), this );
    }

    itkSetMacro( Level, size_t );

protected:
	IterationUpdate(): m_Level(0) {}
	~IterationUpdate(){}

	//void PrintSelf( std::ostream &os, itk::Indent indent ) const;
	size_t m_Level;
private:
	IterationUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented


	OptimizerPointer   m_Optimizer;
};

} // end namespace rstk


#endif /* ITERATIONUPDATE_H_ */
