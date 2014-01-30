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
                            m_UseGridLevelsInitialization(false),
                            m_UseGridSizeInitialization(true),
                            m_Initialized(false) {


}


template < typename TFixedImage, typename TTransform, typename TComputationalValue >
void
ACWERegistrationMethod< TFixedImage, TTransform, TComputationalValue >
::Initialize() {
	if ( m_Initialized ) {
		return;
	}

	ReferenceImagePointer refim = this->GetInput(0);
	GridSizeType maxSize = refim->GetLargestPossibleRegion().GetSize();


	// Planify levels and sizes.
	if ( m_UseGridLevelsInitialization && m_NumberOfLevels>0 ) {
		if ( m_NumberOfLevels == 1 ) {
			m_GridSchedule.push_back( maxSize );
		} else {
			bool invalidRatio = false;

			for( size_t i = 0; i < Dimension; i++) {
				int maxLevels = maxSize[i] - 4;

				if ( maxLevels <= 0 ) {
					itkExceptionMacro(<< "image size must be >= 3 pixels along dimension " << i );
				}

				if ( m_NumberOfLevels > maxLevels ) {
					m_NumberOfLevels = maxLevels;
					itkWarningMacro( << "too many levels required, NumberOfLevels has been updated to " << m_NumberOfLevels );
				}
			}

			m_GridSchedule.resize(m_NumberOfLevels);
			m_GridSchedule[m_NumberOfLevels-1] = maxSize;

			GridSizeType gridStep;
			for( size_t i = 0; i < Dimension; i++){
				gridStep[i] = floor(  1.0*(maxSize[i]-4) / m_NumberOfLevels );
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

	// Initialize LevelSet function
	FunctionalPointer functional = FunctionalType::New();
	functional->SetReferenceImage( im );


	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetFunctional( functional );

	IterationUpdatePointer iup = IterationUpdateType::New();
	iup->SetOptimizer( opt );
	//iup->SetTrackEnergyOn();
#ifndef NDEBUG
	IterationWriterUpdatePointer iwp = IterationWriterUpdate::New();
	iwp->SetOptimizer( opt );
	iwp->SetPrefix( outPrefix );
#endif
}


} // namespace rstk


#endif /* ACWEREGISTRATIONMETHOD_HXX_ */
