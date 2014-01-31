// --------------------------------------------------------------------------------------
// File:          ACWERegistrationMethod.hxx
// Date:          Jan 29, 2014
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
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef ACWEREGISTRATIONMETHOD_HXX_
#define ACWEREGISTRATIONMETHOD_HXX_

#include "ACWERegistrationMethod.h"

#include <boost/lexical_cast.hpp>
#include <algorithm>    // std::fill

namespace rstk {

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::ACWERegistrationMethod(): m_NumberOfLevels(0),
 	 	 	 	 	 	 	m_CurrentLevel(0),
 	 	 	 	 	 	 	m_OutputPrefix(""),
                            m_UseGridLevelsInitialization(false),
                            m_UseGridSizeInitialization(true),
                            m_Initialized(false),
                            m_Stop(false) {
	this->m_StopCondition      = ALL_LEVELS_DONE;
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";

	this->SetNumberOfRequiredOutputs( 1 );

	this->m_OutputTransform = OutputTransformType::New();
	DecoratedOutputTransformPointer transformDecorator = DecoratedOutputTransformType::New().GetPointer();
	transformDecorator->Set( this->m_OutputTransform );
	this->ProcessObject::SetNthOutput( 0, transformDecorator );

	this->m_MinGridSize.Fill( 4 );
}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::PrintSelf( std::ostream &os, itk::Indent indent ) const {
	Superclass::PrintSelf(os,indent);
}


template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::GenerateData() {
	this->InvokeEvent( itk::StartEvent() );

	if ( ! m_Initialized ) {
		this->GenerateSchedule();
	}


	while( this->m_CurrentLevel < this->m_NumberOfLevels ) {
		try {
			this->SetUpLevel( this->m_CurrentLevel );
		} catch ( itk::ExceptionObject & err ) {
			this->Stop( LEVEL_SETTING_ERROR,  "Error while setting up level "
					+  boost::lexical_cast<std::string>(this->m_CurrentLevel) );
			throw err;  // Pass exception to caller
		}

		try {
			m_Optimizers[this->m_CurrentLevel]->Start();
		} catch ( itk::ExceptionObject & err ) {
			this->Stop( LEVEL_PROCESS_ERROR, "Error while executing level "
					+ boost::lexical_cast<std::string>(this->m_CurrentLevel));
			throw err;  // Pass exception to caller
		}
		// Add JSON tree to the general logging facility
		//root["iterations"] = iup->GetJSONRoot();

		this->InvokeEvent( itk::IterationEvent() );
		this->m_CurrentLevel++;

		if ( this->m_CurrentLevel == this->m_NumberOfLevels ) {
			this->Stop( ALL_LEVELS_DONE, "All levels are finished ("
					+ boost::lexical_cast<std::string>(this->m_NumberOfLevels) + " levels)." );
			break;
		}
	}
}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::GenerateSchedule() {
	try {
		ReferenceImageConstPointer refim = this->GetFixedImage();

		for ( size_t i = 0; i<Dimension; i++){
			if ( m_MaxGridSize[i] == 0 ) {
				m_MaxGridSize[i] = refim->GetLargestPossibleRegion().GetSize()[i];
			}

			if ( m_MaxGridSize[i] <= m_MinGridSize[i] ) {
				m_MaxGridSize[i] = m_MinGridSize[i];
			}
		}

		if ( m_MaxGridSize == m_MinGridSize ) {
			this->SetNumberOfLevels( 1 );
		}

		// Schedule levels and sizes.
		if ( m_UseGridLevelsInitialization && m_NumberOfLevels>0 ) {
			if ( m_NumberOfLevels == 1 ) {
				m_GridSchedule.push_back( m_MaxGridSize );
			} else {

				for( size_t i = 0; i < Dimension; i++) {
					int maxLevels = m_MaxGridSize[i] - m_MinGridSize[i];

					if ( maxLevels <= 0 ) {
						itkExceptionMacro(<< "image size must be >= 3 pixels along dimension " << i );
					}

					if ( m_NumberOfLevels > (size_t) maxLevels ) {
						this->SetNumberOfLevels( maxLevels );
						itkWarningMacro( << "too many levels required, NumberOfLevels has been updated to " << m_NumberOfLevels );
					}
				}

				m_GridSchedule[m_NumberOfLevels-1] = m_MaxGridSize;

				GridSizeType gridStep;
				for( size_t i = 0; i < Dimension; i++){
					gridStep[i] = floor(  1.0*(m_MaxGridSize[i]- m_MinGridSize[i]) / (m_NumberOfLevels-1) );
				}

				for( size_t l = m_NumberOfLevels-1; l > 0; --l ){
					GridSizeType prevGrid = m_GridSchedule[l];

					for( size_t i = 0; i < Dimension; i++) {
						prevGrid[i]-= gridStep[i];
					}

					m_GridSchedule[l-1] = prevGrid;
				}
			}
		} // end if m_UseGridLevelsInitialization
		else if ( m_UseGridSizeInitialization ) {

		} else {

		}

		m_Stop = false;
		m_Initialized = true;

	} catch ( itk::ExceptionObject & err ) {
		this->Stop( INITIALIZATION_ERROR, "Error occurred during initialization" );
		throw err;  // Pass exception to caller
	}

	this->InvokeEvent( itk::InitializeEvent() );
}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::SetUpLevel( size_t level ) {
	if( level > (this->m_NumberOfLevels-1) ) {
		itkExceptionMacro( << "Trying to set up a level beyond NumberOfLevels (level=" << (level+1) << ")." );
	}

	ReferenceImageConstPointer im = this->GetFixedImage();

	// Initialize LevelSet function
	m_Functionals[level] = DefaultFunctionalType::New();
	m_Functionals[level]->SetReferenceImage( im );

	if ( level == 0 ) {
		for ( size_t i = 0; i<this->m_Priors.size(); i++ ) {
			m_Functionals[level]->AddShapePrior( this->m_Priors[i] );
		}
	} else {
		for ( size_t i = 0; i<this->m_Priors.size(); i++ ) {
			m_Functionals[level]->AddShapePrior( this->m_Functionals[level-1]->GetCurrentContours()[i] );
		}
	}

	// Connect Optimizer
	m_Optimizers[level] = DefaultOptimizerType::New();
	m_Optimizers[level]->SetFunctional( m_Functionals[level] );

	if ( this->m_NumberOfIterations[level] > 0 ) {
		m_Optimizers[level]->SetNumberOfIterations( this->m_NumberOfIterations[level] );
	}

	// Set registration configuration
	//m_StepSize;
	//m_Alpha;
	//m_Beta;
	//m_DescriptorRecomputationFreq;

//	IterationUpdatePointer iup = IterationUpdateType::New();
//	iup->SetOptimizer( m_Optimizers[level] );
//	//iup->SetTrackEnergyOn();
//#ifndef NDEBUG
//	IterationWriterUpdatePointer iwp = IterationWriterUpdate::New();
//	iwp->SetOptimizer( m_Optimizers[level] );
//	iwp->SetPrefix( m_OutputPrefix );
//#endif

}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::SetNumberOfLevels( size_t levels ) {
	if( levels == 0 || levels > 15 ) {
		itkExceptionMacro( << "intended NumberOfLevels is not valid (zero or >15).")
	}
	this->m_NumberOfLevels = levels;

	m_GridSchedule.resize(m_NumberOfLevels);
	m_Functionals.resize( this->m_NumberOfLevels );
	m_Optimizers.resize( this->m_NumberOfLevels );
	m_NumberOfIterations.resize( this->m_NumberOfLevels );
	m_StepSize.resize( this->m_NumberOfLevels );
	m_Alpha.resize( this->m_NumberOfLevels );
	m_Beta.resize( this->m_NumberOfLevels );
	m_DescriptorRecomputationFreq.resize( this->m_NumberOfLevels );
}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::Stop( StopConditionType code, std::string msg ) {
	this->m_StopConditionDescription << msg;
	this->m_StopCondition = code;
	itkDebugMacro( "Stop called with a description - "  << this->m_StopConditionDescription.str() );
	this->m_Stop = true;
	this->InvokeEvent( itk::EndEvent() );
}

/*
 *  Get output transform
 */
template < typename TFixedImage, typename TTransform, typename TComputationalValue >
const typename ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >::DecoratedOutputTransformType *
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::GetOutput() const
{
  return static_cast<const DecoratedOutputTransformType *>( this->ProcessObject::GetOutput( 0 ) );
}

} // namespace rstk


#endif /* ACWEREGISTRATIONMETHOD_HXX_ */
