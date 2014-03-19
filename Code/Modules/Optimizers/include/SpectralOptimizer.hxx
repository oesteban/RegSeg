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

#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

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
	this->m_Alpha.Fill( 1.0e-2 );
	this->m_Beta.Fill( 1.0e-1 );
	this->m_StopConditionDescription << this->GetNameOfClass() << ": ";
	this->m_GridSize.Fill( 5 );
	this->m_GridSpacing.Fill( 0.0 );
	SplineTransformPointer defaultTransform = SplineTransformType::New();
	this->m_Transform = itkDynamicCastInDebugMode< TransformType* >( defaultTransform.GetPointer() );
	this->m_Transform->SetNumberOfThreads( this->GetNumberOfThreads() );
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
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::ComputeDerivative() {
	// Multiply phi and copy reshaped on this->m_Derivative
	WeightsMatrix phi = this->m_Transform->GetPhi().transpose();
	WeightsMatrix gradVector = this->m_Functional->ComputeDerivative();
	WeightsMatrix derivative( phi.rows(), gradVector.cols() );

	phi.mult( gradVector, derivative );

	typename CoefficientsImageType::PixelType* buff[Dimension];
	for ( size_t i = 0; i<Dimension; i++) {
		this->m_DerivativeCoefficients[i]->FillBuffer( 0.0 );
		buff[i] = this->m_DerivativeCoefficients[i]->GetBufferPointer();
	}
	size_t nPix = this->m_LastField->GetLargestPossibleRegion().GetNumberOfPixels();

	VectorType vi;
	size_t dim;
	VectorType maxSpeed;
	maxSpeed.Fill(0.0);

	typename WeightsMatrix::row row;
	for( size_t r = 0; r<nPix; r++ ){
		vi.Fill(0.0);
		row = derivative.get_row( r );
		for( size_t c = 0; c<row.size(); c++ ) {
			dim = row[c].first;
			vi[dim] = row[c].second;
			*( buff[dim] + r ) = row[c].second;

			if( vi[dim] > maxSpeed[dim] )
				maxSpeed[dim] = vi[dim];
		}
	}

	if( this->m_AutoStepSize && this->m_CurrentIteration < 5 ) {
		this->m_StepSize = (this->m_StepSize  + this->m_LearningRate * ( this->m_MaxDisplacement.GetNorm() / maxSpeed.GetNorm() ) )*0.5;
	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::PostIteration() {
	/* Update the deformation field */
	this->m_Transform->SetCoefficientsImages( this->m_NextCoefficients );
	this->m_Transform->UpdateField();

	this->UpdateField();

	MeasureType meanDisp = this->ComputeIterationChange();

	if ( this->m_UseLightWeightConvergenceChecking ) {
		this->m_CurrentValue = log( meanDisp );
	} else {
		this->m_CurrentValue = this->GetCurrentEnergy();
	}

	this->SetUpdate();

	this->m_Transform->Interpolate();
	this->m_Functional->SetCurrentDisplacements( this->m_Transform->GetOffGridValueMatrix() );
}


template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::GetCurrentRegularizationEnergy() {
	if (!this->m_RegularizationEnergyUpdated ){
		this->m_RegularizationEnergy=0;
		const VectorType* fBuffer = this->m_LastField->GetBufferPointer();
		const VectorType* dBuffer[Dimension];
		size_t nPix = this->m_LastField->GetLargestPossibleRegion().GetNumberOfPixels();

		VectorType u, d_u;
		InternalComputationValueType u_norm = 0.0;
		InternalComputationValueType du_norm = 0.0;


		this->m_Transform->ComputeGradientField();


		double normalizer = 1.0;

		for (size_t i = 0; i<Dimension; i++ ) {
			dBuffer[i] =  this->m_Transform->GetDerivatives()[i]->GetBufferPointer();
			normalizer*= this->m_LastField->GetSpacing()[i];
		}

		for ( size_t pix = 0; pix<nPix; pix++) {
			u = *(fBuffer+pix);
			for ( size_t i = 0; i<Dimension; i++) {
				// Regularization, first term
				this->m_RegularizationEnergy+= this->m_Alpha[i] * u[i] * u[i];

				// Regularization, second term
				d_u = *(dBuffer[i] + pix );
				for( size_t j = 0; j<Dimension; j++) {
					this->m_RegularizationEnergy+= this->m_Beta[j] * d_u[j] * d_u[j];
				}
			}
		}

		this->m_RegularizationEnergy = normalizer * this->m_RegularizationEnergy;
		this->m_RegularizationEnergyUpdated = true;
	}
	return this->m_RegularizationEnergy;
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>::SpectralUpdate(
		CoefficientsImageArray parameters,
		const CoefficientsImageArray lambda,
		CoefficientsImageArray nextParameters,
		bool changeDirection )
{
	itkDebugMacro("Optimizer Spectral Update");

	typename AddFilterType::Pointer add_filter = AddFilterType::New();
	typename MultiplyFilterType::Pointer dir_filter;

	typename MultiplyFilterType::Pointer r_filter = MultiplyFilterType::New();
	r_filter->SetConstant( 1.0/this->m_StepSize );

	//size_t nPix = this->m_Transform->GetNumberOfParameters();

	for( size_t d = 0; d < Dimension; d++ ) {
		r_filter->SetInput( parameters[d] );
		add_filter->SetInput1( r_filter->GetOutput() );

		if (changeDirection) {
			dir_filter = MultiplyFilterType::New();
			dir_filter->SetConstant( -1.0 );
			dir_filter->SetInput( lambda[d] );
			add_filter->SetInput2( dir_filter->GetOutput() );
		} else {
			add_filter->SetInput2( lambda[2] );
		}
		add_filter->Update();

		CoefficientsImagePointer numerator = add_filter->GetOutput();

		FFTPointer fftFilter = FFTType::New();
		fftFilter->SetInput( numerator );
		fftFilter->Update();  // This is required for computing the denominator for first time
		FTDomainPointer fftComponent =  fftFilter->GetOutput();

		this->ApplyRegularizationComponent(d, fftComponent );

		IFFTPointer ifft = IFFTType::New();
		ifft->SetInput( fftComponent );

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
			nextParameters[d],
			ifft->GetOutput()->GetLargestPossibleRegion(),
			nextParameters[d]->GetLargestPossibleRegion()
		);
	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::ApplyRegularizationComponent( size_t d, typename SpectralOptimizer<TFunctional>::FTDomainType* reference ){

	if ( !this->m_DenominatorCached )
		this->InitializeDenominator( reference );

	size_t nPix = this->m_Denominator[d]->GetLargestPossibleRegion().GetNumberOfPixels();

	// Fill reference buffer with elements updated with the filter
	ComplexType* nBuffer = reference->GetBufferPointer();
	PointValueType* dBuffer = this->m_Denominator[d]->GetBufferPointer();

	ComplexType curval;
	ComplexType ddor;
	ComplexType res;

	for (size_t pix = 0; pix < nPix; pix++ ) {
		curval = *(nBuffer+pix);
		ddor = ComplexType( *(dBuffer+pix), 0.0 );
		res = curval * ddor;
		*(nBuffer+pix) = curval * ddor;
	}
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::InitializeDenominator( itk::ImageBase<Dimension>* reference ){
	PointValueType pi2 = 2* vnl_math::pi;
	size_t nPix = reference->GetLargestPossibleRegion().GetNumberOfPixels();
	ControlPointsGridSizeType size = reference->GetLargestPossibleRegion().GetSize();

	for (size_t i = 0; i<Dimension; i++) {
		PointValueType initVal = (1.0/this->m_StepSize) + this->m_Alpha[i];

		CoefficientsImagePointer dimDdor = CoefficientsImageType::New();
		dimDdor->SetSpacing(   reference->GetSpacing() );
		dimDdor->SetDirection( reference->GetDirection() );
		dimDdor->SetRegions(   reference->GetLargestPossibleRegion().GetSize() );
		dimDdor->Allocate();

		PointValueType* buffer = dimDdor->GetBufferPointer();

		if( this->m_Beta.GetNorm() > 0.0 ) {
			// Fill Buffer with denominator data
			PointValueType lag_el; // accumulates the FT{lagrange operator}.
			typename CoefficientsImageType::IndexType idx;
			for (size_t pix = 0; pix < nPix; pix++ ) {
				lag_el = 0.0;
				idx = dimDdor->ComputeIndex( pix );
				for(size_t d = 0; d < Dimension; d++ ) {
					lag_el+= 2.0 * cos( (pi2*idx[d])/size[d] ) - 2.0;
				}
				*(buffer+pix) = 1.0 / (initVal - m_Beta[i] * lag_el);
			}
		}
		else {
			dimDdor->FillBuffer( 1.0 / initVal );
		}


		this->m_Denominator[i] = dimDdor;
	}
	this->m_DenominatorCached = true;
}


template< typename TFunctional >
void
SpectralOptimizer<TFunctional>::UpdateField() {
	itk::ImageAlgorithm::Copy<FieldType, FieldType>(
		this->m_Transform->GetField(), this->m_CurrentDisplacementField,
		this->m_Transform->GetField()->GetLargestPossibleRegion(),
		this->m_CurrentDisplacementField->GetLargestPossibleRegion()
	);


}

template< typename TFunctional >
typename SpectralOptimizer<TFunctional>::MeasureType
SpectralOptimizer<TFunctional>::ComputeIterationChange() {
	const VectorType* fnextBuffer = this->m_CurrentDisplacementField->GetBufferPointer();
	VectorType* fBuffer = this->m_LastField->GetBufferPointer();
	size_t nPix = this->m_LastField->GetLargestPossibleRegion().GetNumberOfPixels();

	InternalComputationValueType totalNorm = 0;
	VectorType t0,t1;
	InternalComputationValueType diff = 0.0;

	this->m_IsDiffeomorphic = true;

	for (size_t pix = 0; pix < nPix; pix++ ) {
		t0 = *(fBuffer+pix);
		t1 = *(fnextBuffer+pix);
		diff = ( t1 - t0 ).GetNorm();
		totalNorm += diff;

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
	return totalNorm/nPix;
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
	/* Initialize next parameters */
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

	this->m_LastField = FieldType::New();
	this->m_LastField->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
	this->m_LastField->SetSpacing(   coeff[0]->GetSpacing() );
	this->m_LastField->SetDirection( coeff[0]->GetDirection() );
	this->m_LastField->SetOrigin(    coeff[0]->GetOrigin() );
	this->m_LastField->Allocate();
	this->m_LastField->FillBuffer( zerov );

	this->m_CurrentDisplacementField = FieldType::New();
	this->m_CurrentDisplacementField->SetRegions(   coeff[0]->GetLargestPossibleRegion() );
	this->m_CurrentDisplacementField->SetSpacing(   coeff[0]->GetSpacing() );
	this->m_CurrentDisplacementField->SetDirection( coeff[0]->GetDirection() );
	this->m_CurrentDisplacementField->SetOrigin(    coeff[0]->GetOrigin() );
	this->m_CurrentDisplacementField->Allocate();
	this->m_CurrentDisplacementField->FillBuffer( zerov );
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
			//("step-size,s", bpo::value< float > (), "step-size value in optimization")
			//("iterations,i", bpo::value< size_t > (), "number of iterations")
			("grid-size,g", bpo::value< size_t > (), "grid size");
			//("update-descriptors,u", bpo::value< size_t > (), "frequency (iterations) to update descriptors of regions (0=no update)");
}

template< typename TFunctional >
void SpectralOptimizer<TFunctional>
::ParseSettings() {
	Superclass::ParseSettings();

	if( this->m_Settings.count( "alpha" ) ) {
		bpo::variable_value v = this->m_Settings["alpha"];
		InternalComputationValueType alpha = v.as<float> ();
		this->SetAlpha( alpha );
	}
	if( this->m_Settings.count( "beta" ) ){
		bpo::variable_value v = this->m_Settings["beta"];
		this->SetBeta( v.as< float >() );
	}
	//if( this->m_Settings.count( "step-size" ) ){
	//	bpo::variable_value v = this->m_Settings["step-size"];
	//	this->SetStepSize( v.as< float >() );
	//}
    //
	//if( this->m_Settings.count( "iterations" ) ){
	//	bpo::variable_value v = this->m_Settings["iterations"];
	//	this->SetNumberOfIterations( v.as< size_t >() );
	//}

	if( this->m_Settings.count( "grid-size" ) ){
		bpo::variable_value v = this->m_Settings["grid-size"];
		this->SetGridSize( v.as< size_t >() );
	}
	//if (this->m_Settings.count("update-descriptors")) {
	//	bpo::variable_value v = this->m_Settings["update-descriptors"];
	//	size_t updDesc =  v.as<size_t>();
	//	this->SetUseDescriptorRecomputation(true);
	//	this->SetDescriptorRecomputationFreq( updDesc );
	//}

}

} // end namespace rstk


#endif /* SPECTRALOPTIMIZER_HXX_ */
