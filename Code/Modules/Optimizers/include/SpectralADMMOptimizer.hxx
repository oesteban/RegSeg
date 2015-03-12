// --------------------------------------------------------------------------------------
// File:          SpectralADMMOptimizer.hxx
// Date:          Jul 31, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
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

#ifndef SPECTRALADMMOPTIMIZER_HXX_
#define SPECTRALADMMOPTIMIZER_HXX_

#include "SpectralADMMOptimizer.h"

using namespace std;

namespace rstk {

/**
 * Default constructor
 */
template< typename TFunctional >
SpectralADMMOptimizer<TFunctional>::SpectralADMMOptimizer() {
	this->m_Rho = 0.80 * (1.0/this->m_StepSize);
}

template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
}


template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>::Iterate() {
	itkDebugMacro("Optimizer Iteration");

	this->UpdateU();
	this->UpdateV();
	this->UpdateLambda();
}

template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>::UpdateU(){
	itkDebugMacro("Optimizer Update u_k");

	ParametersConstPointer shapeGradient = this->m_Functional->GetDerivative();

	VectorType* uFieldBuffer = this->m_NextParameters->GetBufferPointer();
	const VectorType* vFieldBuffer = this->m_vField->GetBufferPointer();
	const VectorType* lFieldBuffer = this->m_lambdaField->GetBufferPointer();
	const VectorType* sFieldBuffer = shapeGradient->GetBufferPointer();
	size_t nPix = this->m_Parameters->GetLargestPossibleRegion().GetNumberOfPixels();

	for ( size_t p = 0; p<nPix; p++ ) {
		*( uFieldBuffer + p ) = *( vFieldBuffer+p ) - *( lFieldBuffer+p ) * this->m_StepSize - *( sFieldBuffer+p ) * this->m_StepSize * this->m_LearningRate;
	}

//#ifndef NDEBUG
//	typedef rstk::DisplacementFieldComponentsFileWriter<DeformationFieldType> Writer;
//	typename Writer::Pointer p = Writer::New();
//	std::stringstream ss;
//	ss << "gradient_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
//	p->SetFileName( ss.str().c_str() );
//	p->SetInput( shapeGradient );
//	p->Update();
//	ss.str("");
//	ss << "nextu_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
//	p->SetFileName( ss.str().c_str() );
//	p->SetInput( this->m_ParametersNext );
//	p->Update();
//#endif

}

template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>::UpdateV() {
	itkDebugMacro("Optimizer Update v");

	this->ComputeUpdate( this->m_NextParameters, this->m_CurrentSpeeds, this->m_vFieldNext, true);

//#ifndef NDEBUG
//	typedef rstk::DisplacementFieldComponentsFileWriter<DeformationFieldType> Writer;
//	typename Writer::Pointer p = Writer::New();
//	std::stringstream ss;
//	ss << "nextv_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
//	p->SetFileName( ss.str().c_str() );
//	p->SetInput( this->m_vFieldNext );
//	p->Update();
//	ss.str("");
//	ss << "nextu2_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
//	p->SetFileName( ss.str().c_str() );
//	p->SetInput( this->m_ParametersNext );
//	p->Update();
//#endif
}

template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>::UpdateLambda() {
	VectorType* lFieldNextBuffer = this->m_lambdaFieldNext->GetBufferPointer();
	const VectorType* uFieldBuffer = this->m_NextParameters->GetBufferPointer();
	const VectorType* vFieldBuffer = this->m_vFieldNext->GetBufferPointer();
	const VectorType* lFieldBuffer = this->m_lambdaField->GetBufferPointer();
	size_t nPix = this->m_NextParameters->GetLargestPossibleRegion().GetNumberOfPixels();

	VectorType diff;
	for ( size_t p = 0; p<nPix; p++ ) {
		diff = ( *( uFieldBuffer+p ) - *( vFieldBuffer+p ) );
		*( lFieldNextBuffer+p ) = *( lFieldBuffer+p ) + diff * this->m_Rho;
	}

//#ifndef NDEBUG
//	typedef rstk::DisplacementFieldComponentsFileWriter<DeformationFieldType> Writer;
//	typename Writer::Pointer p = Writer::New();
//	std::stringstream ss;
//	ss << "nextlambda_" << std::setfill('0')  << std::setw(3) << this->m_CurrentIteration << ".nii.gz";
//	p->SetFileName( ss.str().c_str() );
//	p->SetInput( this->m_lambdaFieldNext );
//	p->Update();
//#endif

}

template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>
::SetUpdate() {
	const VectorType* uNext = this->m_NextParameters->GetBufferPointer();
	const VectorType* vNext = this->m_vFieldNext->GetBufferPointer();
	const VectorType* lNext = this->m_lambdaFieldNext->GetBufferPointer();

	VectorType* uCurr = this->m_Parameters->GetBufferPointer();
	VectorType* vCurr = this->m_vField->GetBufferPointer();
	VectorType* lCurr = this->m_lambdaField->GetBufferPointer();

	size_t nPix = this->m_Parameters->GetLargestPossibleRegion().GetNumberOfPixels();

	for( size_t pix = 0; pix<nPix; pix++){
		*(uCurr+pix) = *(uNext+pix);
		*(vCurr+pix) = *(vNext+pix);
		*(lCurr+pix) = *(lNext+pix);
	}
}


template< typename TFunctional >
void SpectralADMMOptimizer<TFunctional>
::InitializeAuxiliarParameters() {
	VectorType zerov;
	zerov.Fill(0.0);

	typename ParametersType::SpacingType spacing = this->m_Parameters->GetSpacing();
	typename ParametersType::DirectionType dir = this->m_Parameters->GetDirection();
	typename ParametersType::PointType orig = this->m_Parameters->GetOrigin();

	this->m_vField = ParametersType::New();
	this->m_vField->SetRegions( this->m_GridSize );
	this->m_vField->SetSpacing( spacing );
	this->m_vField->SetDirection( dir );
	this->m_vField->SetOrigin( orig );
	this->m_vField->Allocate();
	this->m_vField->FillBuffer( zerov );

	this->m_lambdaField = ParametersType::New();
	this->m_lambdaField->SetRegions( this->m_GridSize );
	this->m_lambdaField->SetSpacing( spacing );
	this->m_lambdaField->SetDirection( dir );
	this->m_lambdaField->SetOrigin( orig );
	this->m_lambdaField->Allocate();
	this->m_lambdaField->FillBuffer( zerov );

	this->m_vFieldNext = ParametersType::New();
	this->m_vFieldNext->SetRegions( this->m_GridSize );
	this->m_vFieldNext->SetSpacing( spacing );
	this->m_vFieldNext->SetDirection( dir );
	this->m_vFieldNext->SetOrigin( orig );
	this->m_vFieldNext->Allocate();
	this->m_vFieldNext->FillBuffer( zerov );

	this->m_lambdaFieldNext = ParametersType::New();
	this->m_lambdaFieldNext->SetRegions( this->m_GridSize );
	this->m_lambdaFieldNext->SetSpacing( spacing );
	this->m_lambdaFieldNext->SetDirection( dir );
	this->m_lambdaFieldNext->SetOrigin( orig );
	this->m_lambdaFieldNext->Allocate();
	this->m_lambdaFieldNext->FillBuffer( zerov );
}

} // end namespace rstk



#endif /* SPECTRALADMMOPTIMIZER_HXX_ */
