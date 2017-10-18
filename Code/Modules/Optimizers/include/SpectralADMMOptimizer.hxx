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
