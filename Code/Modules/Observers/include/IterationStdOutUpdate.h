// --------------------------------------------------------------------------------------
// File:          IterationStdOutUpdate.h
// Date:          Mar 18, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of RegSeg
//
// RegSeg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RegSeg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RegSeg.  If not, see <http://www.gnu.org/licenses/>.
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

    void Execute(const itk::Object * object, const itk::EventObject & event) {
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
