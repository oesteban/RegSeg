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
			this->m_StopCondition = LEVEL_SETTING_ERROR;
			this->m_StopConditionDescription << "Error while setting up level " << this->m_CurrentLevel;
			this->Stop();
			throw err;  // Pass exception to caller
		}

		try {
			m_Optimizers[this->m_CurrentLevel]->Start();
		} catch ( itk::ExceptionObject & err ) {
			this->m_StopCondition = LEVEL_PROCESS_ERROR;
			this->m_StopConditionDescription << "Error while executing level " << this->m_CurrentLevel;
			this->Stop();
			throw err;  // Pass exception to caller
		}
		// Add JSON tree to the general logging facility
		//root["iterations"] = iup->GetJSONRoot();

		this->InvokeEvent( itk::IterationEvent() );
		this->m_CurrentLevel++;

		if ( this->m_CurrentLevel == this->m_NumberOfLevels ) {
			this->m_StopConditionDescription << "All levels are finished (" << this->m_NumberOfLevels << " levels).";
			this->m_StopCondition = ALL_LEVELS_DONE;
			this->Stop();
			break;
		}
	}
}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::GenerateSchedule() {
	try {
		ReferenceImagePointer refim = this->GetInput(0);

		for ( size_t i = 0; i<Dimension; i++){
			if ( m_MaxGridSize[i] == 0 ) {
				m_MaxGridSize[i] = refim->GetLargestPossibleRegion().GetSize()[i];
			}
		}

		// Schedule levels and sizes.
		if ( m_UseGridLevelsInitialization && m_NumberOfLevels>0 ) {
			if ( m_NumberOfLevels == 1 ) {
				m_GridSchedule.push_back( m_MaxGridSize );
			} else {
				bool invalidRatio = false;

				for( size_t i = 0; i < Dimension; i++) {
					int maxLevels = m_MaxGridSize[i] - 4;

					if ( maxLevels <= 0 ) {
						itkExceptionMacro(<< "image size must be >= 3 pixels along dimension " << i );
					}

					if ( m_NumberOfLevels > maxLevels ) {
						m_NumberOfLevels = maxLevels;
						itkWarningMacro( << "too many levels required, NumberOfLevels has been updated to " << m_NumberOfLevels );
					}
				}

				m_GridSchedule.resize(m_NumberOfLevels);
				m_GridSchedule[m_NumberOfLevels-1] = m_MaxGridSize;

				GridSizeType gridStep;
				for( size_t i = 0; i < Dimension; i++){
					gridStep[i] = floor(  1.0*(m_MaxGridSize[i]-4) / m_NumberOfLevels );
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

		m_Functionals.resize( this->m_NumberOfLevels );
		m_Optimizers.resize( this->m_NumberOfLevels );

		m_Stop = false;
		m_Initialized = true;

	} catch ( itk::ExceptionObject & err ) {
		this->m_StopCondition = INITIALIZATION_ERROR;
		this->m_StopConditionDescription << "Error occured during initialization";
		this->Stop();
		throw err;  // Pass exception to caller
	}
}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::SetUpLevel( size_t level ) {
	if( level > (this->m_NumberOfLevels-1) ) {
		itkExceptionMacro( << "Trying to set up a level beyond NumberOfLevels (level=" << (level+1) << ")." );
	}

	ReferenceImageConstPointer im = this->GetInput( 0 );

	// Initialize LevelSet function
	m_Functionals[level] = FunctionalType::New();
	m_Functionals[level]->SetReferenceImage( im );


	// Connect Optimizer
	m_Optimizers[level] = Optimizer::New();
	m_Optimizers[level]->SetFunctional( m_Functionals[level] );

	IterationUpdatePointer iup = IterationUpdateType::New();
	iup->SetOptimizer( m_Optimizers[level] );
	//iup->SetTrackEnergyOn();
#ifndef NDEBUG
	IterationWriterUpdatePointer iwp = IterationWriterUpdate::New();
	iwp->SetOptimizer( m_Optimizers[level] );
	iwp->SetPrefix( m_OutputPrefix );
#endif

}

template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::Stop() {
	itkDebugMacro( "Stop called with a description - "  << this->GetStopConditionDescription() );
	this->m_Stop = true;
	this->InvokeEvent( itk::EndEvent() );
}

} // namespace rstk


#endif /* ACWEREGISTRATIONMETHOD_HXX_ */
