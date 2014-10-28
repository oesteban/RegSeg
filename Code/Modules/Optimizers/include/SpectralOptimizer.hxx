// --------------------------------------------------------------------------------------
// File:          SpectralOptimizer.hxx
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

#ifndef SPECTRALOPTIMIZER_HXX_
#define SPECTRALOPTIMIZER_HXX_

#include "SpectralOptimizer.h"
#include <algorithm>
#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <itkComplexToRealImageFilter.h>
#include <itkImageAlgorithm.h>

using namespace std;


namespace rstk {

/**
 * Default constructor
 */
template< typename TFunctional >
SpectralOptimizer<TFunctional>::SpectralOptimizer():
Superclass(),
m_DenominatorCached( false ),
m_RegularizationEnergy( 0.0 ),
m_CurrentTotalEnergy(itk::NumericTraits<MeasureType>::infinity()),
m_RegularizationEnergyUpdated(true)
{
	this->m_Alpha.Fill( 0.0 );
	this->m_Beta.Fill( 0.0 );
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_GridSize.Fill( 5 );
	this->m_GridSpacing.Fill( 0.0 );
	SplineTransformPointer defaultTransform = SplineTransformType::New();
	this->m_Transform = itkDynamicCastInDebugMode< TransformType* >( defaultTransform.GetPointer() );
	this->m_Transform->SetNumberOfThreads( this->GetNumberOfThreads() );

	SplineTransformPointer initTransform = SplineTransformType::New();
	this->m_TotalTransform = itkDynamicCastInDebugMode< TransformType* >( initTransform.GetPointer() );
	this->m_TotalTransform->SetNumberOfThreads( this->GetNumberOfThreads() );
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
	os << indent << "Learning rate:" << this->m_LearningRate << std::endl;
	os << indent << "Number of iterations: " << this->m_NumberOfIterations << std::endl;
	os << indent << "Current iteration: " << this->m_CurrentIteration << std::endl;
	os << indent << "Stop condition:" << this->m_StopCondition << std::endl;
	os << indent << "Stop condition description: " << this->m_StopConditionDescription.str() << std::endl;
	os << indent << "Alpha: " << this->m_Alpha << std::endl;
	os << indent << "Beta: " << this->m_Beta << std::endl;
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::ComputeDerivative() {
	// Multiply phi and copy reshaped on this->m_Derivative
	WeightsMatrix phi = this->m_Transform->GetPhi().transpose();
	ParametersContainer gradVector = this->m_Functional->ComputeDerivative();
	ParametersContainer derivative;

	typename CoefficientsImageType::PixelType* buff[Dimension];
	for ( size_t i = 0; i<Dimension; i++) {
		phi.mult( gradVector[i], derivative[i] );

		this->m_DerivativeCoefficients[i]->FillBuffer( 0.0 );
		buff[i] = this->m_DerivativeCoefficients[i]->GetBufferPointer();
	}

	size_t nPix = this->m_LastCoeff->GetLargestPossibleRegion().GetNumberOfPixels();
	VectorType vi;
	InternalComputationValueType val;
	size_t dim;
	VectorType maxSpeed, vs;
	maxSpeed.Fill(0.0);

	double norm = this->m_Functional->GetMaxEnergy();
	std::vector< double > speednorms;

	for( size_t r = 0; r<nPix; r++ ){
		vi.Fill(0.0);
		for( size_t c=0; c<Dimension; c++) {
			val = derivative[c][r] / norm;
			vs[c] = val;
			*( buff[c] + r ) = val;
		}
		speednorms.push_back(vs.GetNorm());
	}
	std::sort(speednorms.begin(), speednorms.end());
	this->m_MaximumGradient = speednorms.back();

	//if( this->m_AutoStepSize && this->m_CurrentIteration < 5 ) {
	//	this->m_StepSize = (this->m_StepSize  + this->m_LearningRate * ( this->m_MaxDisplacement.GetNorm() / maxSpeed.GetNorm() ) )*0.5;
	//}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::PostIteration() {
	/* Update the deformation field */
	this->m_Transform->SetCoefficientsImages( this->m_NextCoefficients );
	this->m_Transform->UpdateField();
	//this->UpdateField();

	this->SetUpdate();
	this->ComputeIterationSpeed();

	if ( this->m_UseLightWeightConvergenceChecking ) {
		this->m_CurrentValue = log( 1.0 + this->m_MaxSpeed );
	} else {
		this->m_CurrentValue = this->GetCurrentEnergy();
	}


	this->m_Transform->Interpolate();
	this->m_Functional->SetCurrentDisplacements( this->m_Transform->GetOffGridFieldValues() );
}


template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::GetCurrentRegularizationEnergy() {
	bool haveAlpha = this->m_Alpha.GetNorm() > 1.0e-8;
	bool haveBeta = this->m_Beta.GetNorm() > 1.0e-8;

	if (!haveAlpha && !haveBeta) {
		return 0;
	}

	if (!this->m_RegularizationEnergyUpdated ){
		this->m_RegularizationEnergy=0;

		// Add current coefficients to the corresponding ones at init
		this->m_FieldCoeffAdder->SetInput2(this->m_LastCoeff);
		this->m_FieldCoeffAdder->Update();
		FieldConstPointer currCoeff = this->m_FieldCoeffAdder->GetOutput();

		const VectorType* fBuffer = currCoeff->GetBufferPointer();
		const VectorType* dBuffer;
		FieldConstPointer derivCoeff;

		if (haveBeta) {
			this->m_TotalTransform->SetCoefficientsVectorImage(currCoeff);
			this->m_TotalTransform->ComputeGradientField();
			derivCoeff = this->m_TotalTransform->GetGradientField();
			dBuffer = derivCoeff->GetBufferPointer();
		}

		// Initialize dimensional parameters
		double alpha[Dimension];
		double beta[Dimension];
		double vxvol = 1.0;
		for (size_t i = 0; i<Dimension; i++ ) {
			alpha[i] = this->m_Alpha[i];
			beta[i] = this->m_Beta[i];
			vxvol*= currCoeff->GetSpacing()[i];
		}

		VectorType u;
		VectorType du;
		double u2, du2;
		double energy = 0.0;
		double totalVol = 0.0;
		double eA, eB, e;
		size_t nPix = currCoeff->GetLargestPossibleRegion().GetNumberOfPixels();

		u.Fill(0.0);
		du.Fill(0.0);
		for ( size_t pix = 0; pix<nPix; pix++) {
			if (haveAlpha) u = *(fBuffer+pix);
			if (haveBeta) du = *(dBuffer + pix);

			e = 0.0;
			eA = 0.0;
			eB = 0.0;
			for ( size_t i = 0; i<Dimension; i++) {
				// Regularization, first term
				u2 = u[i] * u[i];
				if (u2 > 1.0e-4)
					eA+= alpha[i] * u2;

				// Regularization, second term
				du2 = du[i] * du[i];
				if (du2 > 1.0e-5)
					eB+= beta[i] * du2;
			}
			energy+= (eA + eB);
		}

		this->m_RegularizationEnergy = energy;
		this->m_RegularizationEnergyUpdated = true;
	}
	return this->m_RegularizationEnergy;
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::ComputeUpdate(
		CoefficientsImageArray uk,
		const CoefficientsImageArray gk,
		CoefficientsImageArray next_uk,
		bool changeDirection){

	typename AddFilterType::Pointer add_f = AddFilterType::New();
	typename MultiplyFilterType::Pointer dir_filter;

	typename MultiplyFilterType::Pointer step_f = MultiplyFilterType::New();
	if (changeDirection) {
		step_f->SetConstant( -1.0 * this->m_StepSize );
	} else {
		step_f->SetConstant( this->m_StepSize );
	}

	typename MultiplyFilterType::Pointer scale_f = MultiplyFilterType::New();
	CoefficientsImagePointer result;
	InternalVectorType s;

	for( size_t d = 0; d < Dimension; d++ ) {
		step_f->SetInput( gk[d] );
		step_f->Update();
		add_f->SetInput1( uk[d] );
		add_f->SetInput2( step_f->GetOutput() );
		add_f->Update();
		s[d] = 1.0;

		if( this->m_Alpha[d] > 1.0e-8) {
			s[d] = 1.0 / (1.0 + 2.0 * this->m_Alpha[d] * this->m_StepSize);
			scale_f->SetConstant( s[d] );
			scale_f->SetInput(add_f->GetOutput());
			scale_f->Update();
			result = scale_f->GetOutput();
		} else {
			result = add_f->GetOutput();
		}

		if( this->m_Beta[d] > 1.0e-8) {
			InternalComputationValueType scaler = 2.0 * this->m_Beta[d] * this->m_StepSize * s[d];
			this->BetaRegularization(result, next_uk, scaler, d);
		} else {
			// Set component in destination buffer of coefficients
			itk::ImageAlgorithm::Copy< CoefficientsImageType, CoefficientsImageType >(
				result, next_uk[d],
				result->GetLargestPossibleRegion(),
				next_uk[d]->GetLargestPossibleRegion()
			);
		}
	}

//	double norm = 0.0;
//	double maxs = 0.0;
//	VectorType v;
//	const typename CoefficientsImageType::PixelType* buffer[3];
//	size_t nPix = next_uk[0]->GetLargestPossibleRegion().GetNumberOfPixels();
//	for( size_t d = 0; d < Dimension; d++ ) {
//		buffer[d] = next_uk[d]->GetBufferPointer();
//	}
//	for(size_t i = 0; i < nPix; i++) {
//		for(size_t d = 0; d < Dimension; d++) {
//			v[d] = *(buffer[d] + i);
//		}
//		norm = v.GetNorm();
//		if (norm > maxs) {
//			maxs = norm;
//		}
//	}
//
//	std::cout << "Speed=" << maxs << "." << std::endl;
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::BetaRegularization(
		CoefficientsImagePointer numerator,
		CoefficientsImageArray next_uk,
		InternalComputationValueType s,
		size_t d) {
	itkDebugMacro("Optimizer Spectral Update");

	FFTPointer fftFilter = FFTType::New();
	fftFilter->SetInput( numerator );
	fftFilter->Update();  // This is required for computing the denominator for first time

	if ( !this->m_DenominatorCached )
		this->InitializeDenominator( fftFilter->GetOutput() );

	typename FTMultiplyFilterType::Pointer scale_f = FTMultiplyFilterType::New();
	scale_f->SetInput(this->m_FTLaplacian);
	scale_f->SetConstant( -1.0 * s);

	typename FTAddFilterType::Pointer add_f = FTAddFilterType::New();
	add_f->SetInput1(this->m_FTOnes);
	add_f->SetInput2(scale_f->GetOutput());

	typename FTDivideFilterType::Pointer div_f = FTDivideFilterType::New();
	div_f->SetInput1(fftFilter->GetOutput());
	div_f->SetInput2(add_f->GetOutput());
	div_f->Update();

	IFFTPointer ifft = IFFTType::New();
	ifft->SetInput( div_f->GetOutput() );

	try {
		ifft->Update();
	}
	catch ( itk::ExceptionObject & err ) {
		this->m_StopCondition = Superclass::UPDATE_PARAMETERS_ERROR;
		this->m_StopConditionDescription << "Optimizer error: spectral update failed";
		this->Stop();
		// Pass exception to caller
		throw err;
	}

	// Set component in destination buffer of coefficients
	itk::ImageAlgorithm::Copy< CoefficientsImageType, CoefficientsImageType >(
		ifft->GetOutput(),
		next_uk[d],
		ifft->GetOutput()->GetLargestPossibleRegion(),
		next_uk[d]->GetLargestPossibleRegion()
	);

	// typedef itk::ComplexToRealImageFilter<FTDomainType, CoefficientsImageType> RealFilterType;
	// typename RealFilterType::Pointer realFilter = RealFilterType::New();
	// realFilter->SetInput(add_f->GetOutput());
	// realFilter->Update();
    //
	// std::stringstream ss;
	// ss << "test_coeff" << d << "_it" << this->m_CurrentIteration << ".nii.gz";
	// typename itk::ImageFileWriter<CoefficientsImageType>::Pointer w = itk::ImageFileWriter<CoefficientsImageType>::New();
	// w->SetFileName(ss.str().c_str());
	// w->SetInput(realFilter->GetOutput());
	// w->Update();
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::InitializeDenominator( itk::ImageBase<Dimension>* reference ){
	PointValueType pi2 = 2* vnl_math::pi;
	size_t nPix = reference->GetLargestPossibleRegion().GetNumberOfPixels();
	ControlPointsGridSizeType size = reference->GetLargestPossibleRegion().GetSize();

	ComplexType one = itk::NumericTraits<ComplexType>::One;
	this->m_FTOnes = FTDomainType::New();
	this->m_FTOnes->SetSpacing(   reference->GetSpacing() );
	this->m_FTOnes->SetDirection( reference->GetDirection() );
	this->m_FTOnes->SetRegions(   reference->GetLargestPossibleRegion().GetSize() );
	this->m_FTOnes->Allocate();
	this->m_FTOnes->FillBuffer(one);

	this->m_FTLaplacian = FTDomainType::New();
	this->m_FTLaplacian->SetSpacing(   reference->GetSpacing() );
	this->m_FTLaplacian->SetDirection( reference->GetDirection() );
	this->m_FTLaplacian->SetRegions(   reference->GetLargestPossibleRegion().GetSize() );
	this->m_FTLaplacian->Allocate();
	this->m_FTLaplacian->FillBuffer( itk::NumericTraits<ComplexType>::Zero);

	// Fill Buffer with denominator data
	ComplexType* buffer = this->m_FTLaplacian->GetBufferPointer();
	InternalComputationValueType lag_el; // accumulates the FT{lagrange operator}.
	typename FTDomainType::IndexType idx;
	for (size_t pix = 0; pix < nPix; pix++ ) {
		lag_el = 0.0;
		idx = this->m_FTLaplacian->ComputeIndex(pix);
		for(size_t d = 0; d < Dimension; d++ ) {
			lag_el+= 2.0 * cos( (pi2*idx[d])/size[d] ) - 2.0;
		}
		*(buffer+pix) = one * lag_el;
	}
	this->m_DenominatorCached = true;
}


template< typename TFunctional >
void
SpectralOptimizer<TFunctional>::UpdateField() {
	itk::ImageAlgorithm::Copy<FieldType, FieldType>(
		this->m_Transform->GetField(), this->m_CurrentCoefficients,
		this->m_Transform->GetField()->GetLargestPossibleRegion(),
		this->m_CurrentCoefficients->GetLargestPossibleRegion()
	);


}

template< typename TFunctional >
void
SpectralOptimizer<TFunctional>::ComputeIterationSpeed() {
	const VectorType* fnextBuffer = this->m_CurrentCoefficients->GetBufferPointer();
	VectorType* fBuffer = this->m_LastCoeff->GetBufferPointer();
	size_t nPix = this->m_LastCoeff->GetLargestPossibleRegion().GetNumberOfPixels();

	InternalComputationValueType totalNorm = 0;
	VectorType t0,t1;
	InternalComputationValueType diff = 0.0;

	this->m_IsDiffeomorphic = true;
	std::vector< InternalComputationValueType > speednorms;
	for (size_t pix = 0; pix < nPix; pix++ ) {
		t0 = *(fBuffer+pix);
		t1 = *(fnextBuffer+pix);
		diff = ( t1 - t0 ).GetNorm();
		totalNorm += diff;
		speednorms.push_back(diff);
		*(fBuffer+pix) = t1; // Copy current to last, once evaluated

		if ( this->m_IsDiffeomorphic ) {
			for( size_t i = 0; i<Dimension; i++) {
				if ( fabs(t1[i]) > this->m_MaxDisplacement[i] ) {
					this->m_IsDiffeomorphic = false;
				}
			}
		}
	}

	this->m_RegularizationEnergyUpdated = (totalNorm==0);

	std::sort(speednorms.begin(), speednorms.end());
	this->m_MaxSpeed = speednorms.back();
	this->m_MeanSpeed = speednorms[int(0.5*(speednorms.size()-1))];
	this->m_AvgSpeed = totalNorm / nPix;
}


template< typename TFunctional >
void SpectralOptimizer<TFunctional>::InitializeParameters() {
	// Check functional exists and hold a reference image
	if ( this->m_Functional.IsNull() ) {
		itkExceptionMacro( << "functional must be set." );
	}

	this->m_Transform->SetControlPointsSize( this->m_GridSize );
	this->m_Transform->SetPhysicalDomainInformation( this->m_Functional->GetReferenceImage() );
	this->m_Transform->SetOffGridPositions( this->m_Functional->GetNodesPosition() );

	CoefficientsImageArray coeff = this->m_Transform->GetCoefficientsImages();
	this->m_GridSpacing = coeff[0]->GetSpacing();

	for (size_t i = 0; i<Dimension; i++) {
		this->m_MaxDisplacement[i] = 0.40 * this->m_GridSpacing[i];
	}


	VectorType zerov; zerov.Fill( 0.0 );
	/* Initialize next uk */
	for ( size_t i=0; i<coeff.Size(); i++ ) {
		this->m_Coefficients[i] = CoefficientsImageType::New();
		this->m_Coefficients[i]->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
		this->m_Coefficients[i]->SetSpacing(   coeff[0]->GetSpacing() );
		this->m_Coefficients[i]->SetDirection( coeff[0]->GetDirection() );
		this->m_Coefficients[i]->SetOrigin(    coeff[0]->GetOrigin() );
		this->m_Coefficients[i]->Allocate();
		this->m_Coefficients[i]->FillBuffer( 0.0 );

		this->m_NextCoefficients[i] = CoefficientsImageType::New();
		this->m_NextCoefficients[i]->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
		this->m_NextCoefficients[i]->SetSpacing(   coeff[0]->GetSpacing() );
		this->m_NextCoefficients[i]->SetDirection( coeff[0]->GetDirection() );
		this->m_NextCoefficients[i]->SetOrigin(    coeff[0]->GetOrigin() );
		this->m_NextCoefficients[i]->Allocate();
		this->m_NextCoefficients[i]->FillBuffer( 0.0 );

		this->m_DerivativeCoefficients[i] = CoefficientsImageType::New();
		this->m_DerivativeCoefficients[i]->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
		this->m_DerivativeCoefficients[i]->SetSpacing(   coeff[0]->GetSpacing() );
		this->m_DerivativeCoefficients[i]->SetDirection( coeff[0]->GetDirection() );
		this->m_DerivativeCoefficients[i]->SetOrigin(    coeff[0]->GetOrigin() );
		this->m_DerivativeCoefficients[i]->Allocate();
		this->m_DerivativeCoefficients[i]->FillBuffer( 0.0 );
	}

	this->m_LastCoeff = FieldType::New();
	this->m_LastCoeff->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
	this->m_LastCoeff->SetSpacing(   coeff[0]->GetSpacing() );
	this->m_LastCoeff->SetDirection( coeff[0]->GetDirection() );
	this->m_LastCoeff->SetOrigin(    coeff[0]->GetOrigin() );
	this->m_LastCoeff->Allocate();
	this->m_LastCoeff->FillBuffer( zerov );

	this->m_CurrentCoefficients = FieldType::New();
	this->m_CurrentCoefficients->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
	this->m_CurrentCoefficients->SetSpacing(   coeff[0]->GetSpacing() );
	this->m_CurrentCoefficients->SetDirection( coeff[0]->GetDirection() );
	this->m_CurrentCoefficients->SetOrigin(    coeff[0]->GetOrigin() );
	this->m_CurrentCoefficients->Allocate();
	this->m_CurrentCoefficients->FillBuffer( zerov );

	if( this->m_InitialDisplacementField.IsNull() ) {
		this->m_InitialDisplacementField = FieldType::New();
		this->m_InitialDisplacementField->SetRegions(   this->m_Functional->GetReferenceImage()->GetLargestPossibleRegion() );
		this->m_InitialDisplacementField->SetSpacing(   this->m_Functional->GetReferenceImage()->GetSpacing() );
		this->m_InitialDisplacementField->SetDirection( this->m_Functional->GetReferenceImage()->GetDirection() );
		this->m_InitialDisplacementField->SetOrigin(    this->m_Functional->GetReferenceImage()->GetOrigin() );
		this->m_InitialDisplacementField->Allocate();
		this->m_InitialDisplacementField->FillBuffer( zerov );
	}

	this->m_TotalTransform->CopyGridInformation( coeff[0] );
	this->m_TotalTransform->SetOutputReference(this->m_InitialDisplacementField);
	this->m_TotalTransform->SetField(this->m_InitialDisplacementField);
	this->m_TotalTransform->ComputeCoefficients();

	this->m_InitialCoeff = FieldType::New();
	this->m_InitialCoeff->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
	this->m_InitialCoeff->SetSpacing(   coeff[0]->GetSpacing() );
	this->m_InitialCoeff->SetDirection( coeff[0]->GetDirection() );
	this->m_InitialCoeff->SetOrigin(    coeff[0]->GetOrigin() );
	this->m_InitialCoeff->Allocate();
	this->m_InitialCoeff->FillBuffer( zerov );

	itk::ImageAlgorithm::Copy<FieldType, FieldType>(this->m_TotalTransform->GetCoefficientsVectorImage(),
			this->m_InitialCoeff, this->m_TotalTransform->GetCoefficientsVectorImage()->GetLargestPossibleRegion(),
			this->m_InitialCoeff->GetLargestPossibleRegion());

	this->m_FieldCoeffAdder = AddFieldFilterType::New();
	this->m_FieldCoeffAdder->SetInput1(this->m_InitialCoeff);
}


template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::GetCurrentEnergy() {
	this->m_CurrentTotalEnergy = this->m_Functional->GetValue() + this->GetCurrentRegularizationEnergy();
	return this->m_CurrentTotalEnergy;
}


template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::AddOptions( SettingsDesc& opts ) {
	Superclass::AddOptions( opts );
	opts.add_options()
			("alpha,a", bpo::value< float > (), "alpha value in regularization")
			("beta,b", bpo::value< float > (), "beta value in regularization")
			("grid-size,g", bpo::value< size_t > (), "grid size");
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::ParseSettings() {
	Superclass::ParseSettings();

	if( this->m_Settings.count( "alpha" ) ) {
		bpo::variable_value v = this->m_Settings["alpha"];
		std::vector<float> alpha = v.as< std::vector<float> > ();
		if (alpha.size() == 1) {
			this->SetAlpha( alpha[0] );
		} else if (alpha.size() == Dimension) {
			InternalVectorType v;
			for( size_t i = 0; i < Dimension; i++)
				v[i] = alpha[i];
			this->SetAlpha(v);
		}
	}
	if( this->m_Settings.count( "beta" ) ){
		bpo::variable_value v = this->m_Settings["beta"];
		std::vector<float> beta = v.as< std::vector<float> > ();
		if (beta.size() == 1) {
			this->SetBeta( beta[0] );
		} else if (beta.size() == Dimension) {
			InternalVectorType v;
			for( size_t i = 0; i < Dimension; i++)
				v[i] = beta[i];
			this->SetBeta(v);
		}
	}

	if( this->m_Settings.count( "grid-size" ) ){
		bpo::variable_value v = this->m_Settings["grid-size"];
		this->SetGridSize( v.as< size_t >() );
	}
}

} // end namespace rstk


#endif /* SPECTRALOPTIMIZER_HXX_ */
