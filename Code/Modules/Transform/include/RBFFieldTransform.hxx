/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __RBFFieldTransform_hxx
#define __RBFFieldTransform_hxx

#include "RBFFieldTransform.h"
#include <itkVectorLinearInterpolateImageFunction.h>

#include <itkImageRegionIteratorWithIndex.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>

namespace rstk {

/**
 * Constructor
 */
template<typename TScalar, unsigned int NDimensions>
RBFFieldTransform<TScalar, NDimensions>::RBFFieldTransform() :
		Superclass(0) {
	this->m_DisplacementField = ITK_NULLPTR;
	this->m_CoefficientsField = ITK_NULLPTR;
	this->m_InverseDisplacementField = ITK_NULLPTR;

	SizeValueType elems = NDimensions * (NDimensions + 3);
	this->m_FixedParameters.SetSize(elems * 2);
	this->m_FixedParameters.Fill(0.0);

	// Setup and assign default interpolator
	typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType,
			ScalarType> DefaultInterpolatorType;
	typename DefaultInterpolatorType::Pointer interpolator =
			DefaultInterpolatorType::New();
	this->m_Interpolator = interpolator;

	typename DefaultInterpolatorType::Pointer inverseInterpolator =
			DefaultInterpolatorType::New();
	this->m_InverseInterpolator = inverseInterpolator;

	// Setup and assign parameter helper. This will hold the displacement field
	// for access through the common OptimizerParameters interface.
	OptimizerParametersHelperType* helper = new OptimizerParametersHelperType;
	// After assigning this, m_Parametes will manage this,
	// deleting when appropriate.
	this->m_Parameters.SetHelper(helper);

	m_DisplacementFieldSetTime = 0;

	/* Initialize the identity jacobian. */
	m_IdentityJacobian.SetSize(NDimensions, NDimensions);
	m_IdentityJacobian.Fill(0.0);
	for (unsigned int dim = 0; dim < NDimensions; dim++) {
		m_IdentityJacobian[dim][dim] = 1.0;
	}
}

/**
 * Destructor
 */
template<typename TScalar, unsigned int NDimensions>
RBFFieldTransform<TScalar, NDimensions>::~RBFFieldTransform() {
}

/**
 * Transform point
 */
template<typename TScalar, unsigned int NDimensions>
typename RBFFieldTransform<TScalar, NDimensions>::OutputPointType RBFFieldTransform<
		TScalar, NDimensions>::TransformPoint(
		const InputPointType& inputPoint) const {
	if (!this->m_DisplacementField) {
		itkExceptionMacro("No displacement field is specified.");
	}
	if (!this->m_CoefficientsField) {
		itkExceptionMacro("No coefficients field is specified.");
	}
	if (!this->m_Interpolator) {
		itkExceptionMacro("No interpolator is specified.");
	}

	typename InterpolatorType::ContinuousIndexType cidx;
	typename InterpolatorType::PointType point;
	point.CastFrom(inputPoint);

	OutputPointType outputPoint;
	outputPoint.CastFrom(inputPoint);

	if (this->m_Interpolator->IsInsideBuffer(point)) {
		this->m_DisplacementField->TransformPhysicalPointToContinuousIndex(
				point, cidx);
		typename InterpolatorType::OutputType displacement =
				this->m_Interpolator->EvaluateAtContinuousIndex(cidx);
		for (unsigned int ii = 0; ii < NDimensions; ++ii) {
			outputPoint[ii] += displacement[ii];
		}
	}
	// else
	// simply return inputPoint

	return outputPoint;
}

/**
 * return an inverse transformation
 */
template<typename TScalar, unsigned int NDimensions>
bool RBFFieldTransform<TScalar, NDimensions>::GetInverse(
		Self *inverse) const {
	if (!inverse || !this->m_InverseDisplacementField) {
		return false;
	} else {
		inverse->SetDisplacementField(this->m_InverseDisplacementField);
		inverse->SetInverseDisplacementField(this->m_DisplacementField);
		inverse->SetInterpolator(this->m_InverseInterpolator);
		inverse->SetInverseInterpolator(this->m_Interpolator);

		return true;
	}
}

// Return an inverse of this transform
template<typename TScalar, unsigned int NDimensions>
typename RBFFieldTransform<TScalar, NDimensions>::InverseTransformBasePointer RBFFieldTransform<
		TScalar, NDimensions>::GetInverseTransform() const {
	Pointer inverseTransform = New();

	if (this->GetInverse(inverseTransform)) {
		return inverseTransform.GetPointer();
	} else {
		return ITK_NULLPTR;
	}
}

/*
 * ComputeJacobianWithRespectToParameters methods
 */

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::ComputeJacobianWithRespectToPosition(
		const InputPointType & point, JacobianType & jacobian) const {
	IndexType idx;

	this->m_DisplacementField->TransformPhysicalPointToIndex(point, idx);
	this->ComputeJacobianWithRespectToPosition(idx, jacobian);
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::ComputeJacobianWithRespectToPosition(
		const IndexType & index, JacobianType & jacobian) const {
	this->ComputeJacobianWithRespectToPositionInternal(index, jacobian, false);
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::ComputeInverseJacobianWithRespectToPosition(
		const InputPointType & point, JacobianType & jacobian) const {
	IndexType idx;
	this->m_DisplacementField->TransformPhysicalPointToIndex(point, idx);
	this->ComputeJacobianWithRespectToPositionInternal(idx, jacobian, true);
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::GetInverseJacobianOfForwardFieldWithRespectToPosition(
		const InputPointType & point, JacobianType & jacobian,
		bool useSVD) const {
	IndexType idx;

	this->m_DisplacementField->TransformPhysicalPointToIndex(point, idx);
	this->GetInverseJacobianOfForwardFieldWithRespectToPosition(idx, jacobian,
			useSVD);
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::GetInverseJacobianOfForwardFieldWithRespectToPosition(
		const IndexType & index, JacobianType & jacobian, bool useSVD) const {
	if (useSVD) {
		this->ComputeJacobianWithRespectToPositionInternal(index, jacobian,
				false);
		vnl_svd<typename JacobianType::ValueType> svd(jacobian);
		for (unsigned int i = 0; i < jacobian.rows(); i++) {
			for (unsigned int j = 0; j < jacobian.cols(); j++) {
				jacobian(i, j) = svd.inverse()(i, j);
			}
		}
	} else {
		this->ComputeJacobianWithRespectToPositionInternal(index, jacobian,
				true);
	}
}

/*
 * ComputeJacobianWithRespectToPositionInternal. Worker method.
 */
template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::ComputeJacobianWithRespectToPositionInternal(
		const IndexType & index, JacobianType & jacobian,
		bool doInverseJacobian) const {
	jacobian.SetSize(NDimensions, NDimensions);

	typename DisplacementFieldType::SizeType size =
			this->m_DisplacementField->GetLargestPossibleRegion().GetSize();
	typename DisplacementFieldType::SpacingType spacing =
			this->m_DisplacementField->GetSpacing();

	IndexType ddrindex;
	IndexType ddlindex;
	IndexType difIndex[NDimensions][2];

	// index offset
	unsigned int posoff = itk::NumericTraits<unsigned int>::One;

	// space between indices
	TScalar space = itk::NumericTraits<TScalar>::One;

	// minimum distance between neighbors
	TScalar mindist = itk::NumericTraits<TScalar>::One;

	// flag indicating a valid location for jacobian calculation
	bool oktosample = true;

	// multiplier for getting inverse jacobian
	TScalar dPixSign = itk::NumericTraits<TScalar>::One;
	dPixSign = doInverseJacobian ? -dPixSign : dPixSign;
	for (unsigned int row = 0; row < NDimensions; row++) {
		TScalar dist = fabs((float) index[row]);
		if (dist < mindist) {
			oktosample = false;
		}
		dist = fabs((TScalar) size[row] - (TScalar) index[row]);
		if (dist < mindist) {
			oktosample = false;
		}
	}

	if (oktosample) {
		OutputVectorType cpix = this->m_DisplacementField->GetPixel(index);
		m_DisplacementField->TransformLocalVectorToPhysicalVector(cpix, cpix);
		// cpix = directionRaw->TransformVector( cpix );
		// itkCentralDifferenceImageFunction does not support 4th order so
		// do manually here
		for (unsigned int row = 0; row < NDimensions; row++) {
			difIndex[row][0] = index;
			difIndex[row][1] = index;
			ddrindex = index;
			ddlindex = index;
			if ((int) index[row] < (int) (size[row] - 2)) {
				difIndex[row][0][row] = index[row] + posoff;
				ddrindex[row] = index[row] + posoff * 2;
			}
			if (index[row] > 1) {
				difIndex[row][1][row] = index[row] - 1;
				ddlindex[row] = index[row] - 2;
			}

			OutputVectorType rpix = m_DisplacementField->GetPixel(
					difIndex[row][1]);
			OutputVectorType lpix = m_DisplacementField->GetPixel(
					difIndex[row][0]);
			OutputVectorType rrpix = m_DisplacementField->GetPixel(ddrindex);
			OutputVectorType llpix = m_DisplacementField->GetPixel(ddlindex);

			m_DisplacementField->TransformLocalVectorToPhysicalVector(rpix,
					rpix);
			m_DisplacementField->TransformLocalVectorToPhysicalVector(rrpix,
					rrpix);
			m_DisplacementField->TransformLocalVectorToPhysicalVector(lpix,
					lpix);
			m_DisplacementField->TransformLocalVectorToPhysicalVector(llpix,
					llpix);

			// 4th order centered difference
			OutputVectorType dPix = (lpix * 8.0 + llpix - rrpix - rpix * 8.0)
					* space / (12.0) * dPixSign;

			// typename DisplacementFieldType::PixelType dPix=
			//      ( lpix - rpix )*space/(2.0*h); //2nd order centered difference
			for (unsigned int col = 0; col < NDimensions; col++) {
				TScalar val = dPix[col] / spacing[col];
				if (row == col) {
					val += 1.0;
				}
				jacobian(col, row) = val;
				// Verify it's a real number
				if (!vnl_math_isfinite(val)) {
					oktosample = false;
					break;
				}
			}
		} // for row
	}   // if oktosample

	if (!oktosample) {
		jacobian.Fill(0.0);
		for (unsigned int i = 0; i < NDimensions; i++) {
			jacobian(i, i) = 1.0;
		}
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::UpdateTransformParameters(
		const DerivativeType & update, ScalarType factor) {
	// This simply adds the values.
	// TODO: This should be multi-threaded probably, via image add filter.
	Superclass::UpdateTransformParameters(update, factor);
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetDisplacementField(
		DisplacementFieldType* field) {
	itkDebugMacro("setting DisplacementField to " << field);
	if (this->m_DisplacementField != field) {
		this->m_DisplacementField = field;

		if (!this->m_InverseDisplacementField.IsNull()) {
			this->m_InverseDisplacementField = ITK_NULLPTR;
		}
		this->Modified();
		/* Store this separately for use in smoothing because we only want
		 * to know when the displacement field object has changed, not just
		 * its contents. */
		this->m_DisplacementFieldSetTime = this->GetMTime();
		if (!this->m_Interpolator.IsNull()) {
			this->m_Interpolator->SetInputImage(this->m_DisplacementField);
		}
	}
	this->SetFixedParametersFromDisplacementField();
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetCoefficientsField(
		DisplacementFieldType* field) {
	itkDebugMacro("setting DisplacementField to " << field);
	if (this->m_CoefficientsField != field) {
		this->m_CoefficientsField = field;
		this->Modified();
		/* Store this separately for use in smoothing because we only want
		 * to know when the displacement field object has changed, not just
		 * its contents. */
		this->m_CoefficientsFieldSetTime = this->GetMTime();
		if (!this->m_Interpolator.IsNull()) {
			this->m_Interpolator->SetInputImage(this->m_CoefficientsField);
		}
		// Assign to parameters object
		this->m_Parameters.SetParametersObject(this->m_CoefficientsField);
	}
	this->SetFixedParametersFromCoefficientsField();
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetInverseDisplacementField(
		DisplacementFieldType* inverseField) {
	itkDebugMacro("setting InverseDisplacementField to " << inverseField);
	if (this->m_InverseDisplacementField != inverseField) {
		this->m_InverseDisplacementField = inverseField;
		if (!this->m_DisplacementField.IsNull() && inverseField) {
			this->VerifyFixedParametersInformation();
		}
		if (!this->m_InverseInterpolator.IsNull()) {
			this->m_InverseInterpolator->SetInputImage(
					this->m_InverseDisplacementField);
		}
		this->Modified();
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::VerifyFixedParametersInformation() {
	if (!this->m_DisplacementField.IsNull()
			&& !this->m_InverseDisplacementField.IsNull()) {
		// check to see if the candidate inverse displacement field has the
		// same fixed parameters as the displacement field.

		SizeType inverseFieldSize =
				this->m_InverseDisplacementField->GetLargestPossibleRegion().GetSize();
		PointType inverseFieldOrigin =
				this->m_InverseDisplacementField->GetOrigin();
		SpacingType inverseFieldSpacing =
				this->m_InverseDisplacementField->GetSpacing();
		DirectionType inverseFieldDirection =
				this->m_InverseDisplacementField->GetDirection();

		SizeType fieldSize =
				this->m_DisplacementField->GetLargestPossibleRegion().GetSize();
		PointType fieldOrigin = this->m_DisplacementField->GetOrigin();
		SpacingType fieldSpacing = this->m_DisplacementField->GetSpacing();
		DirectionType fieldDirection =
				this->m_DisplacementField->GetDirection();

		// tolerance for origin and spacing depends on the size of pixel
		// tolerance for directions a fraction of the unit cube.
		const double coordinateTolerance = 1.0e-6 * fieldSpacing[0];
		const double directionTolerance = 1.0e-6;

		std::ostringstream sizeString;
		std::ostringstream originString;
		std::ostringstream spacingString;
		std::ostringstream directionString;

		bool unequalSizes = false;
		bool unequalOrigins = false;
		bool unequalSpacings = false;
		bool unequalDirections = false;

		if (inverseFieldSize != fieldSize) {
			unequalSizes = true;
			sizeString << "InverseDisplacementField Size: " << inverseFieldSize
					<< ", DisplacementField Size: " << fieldSize << std::endl;

		}
		if (!inverseFieldOrigin.GetVnlVector().is_equal(
				fieldOrigin.GetVnlVector(), coordinateTolerance)) {
			unequalOrigins = true;
			originString << "InverseDisplacementField Origin: "
					<< inverseFieldOrigin << ", DisplacementField Origin: "
					<< fieldOrigin << std::endl;

		}
		if (!inverseFieldSpacing.GetVnlVector().is_equal(
				fieldSpacing.GetVnlVector(), coordinateTolerance)) {
			unequalSpacings = false;
			originString << "InverseDisplacementField Spacing: "
					<< inverseFieldSpacing << ", DisplacementField Spacing: "
					<< fieldSpacing << std::endl;

		}
		if (!inverseFieldDirection.GetVnlMatrix().as_ref().is_equal(
				fieldDirection.GetVnlMatrix(), directionTolerance)) {
			unequalDirections = true;
			originString << "InverseDisplacementField Direction: "
					<< inverseFieldDirection
					<< ", DisplacementField Direction: " << fieldDirection
					<< std::endl;

		}
		if (unequalSizes || unequalOrigins || unequalSpacings
				|| unequalDirections) {
			itkExceptionMacro(
					"The inverse and displacement fields do not have the same fixed parameters: " << std::endl << sizeString.str() << originString.str() << spacingString.str() << directionString.str());
		}
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetInterpolator(
		InterpolatorType* interpolator) {
	itkDebugMacro("setting Interpolator to " << interpolator);
	if (this->m_Interpolator != interpolator) {
		this->m_Interpolator = interpolator;
		this->Modified();
		if (!this->m_DisplacementField.IsNull()) {
			this->m_Interpolator->SetInputImage(this->m_DisplacementField);
		}
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetInverseInterpolator(
		InterpolatorType* interpolator) {
	itkDebugMacro("setting InverseInterpolator to " << interpolator);
	if (this->m_InverseInterpolator != interpolator) {
		this->m_InverseInterpolator = interpolator;
		this->Modified();
		if (!this->m_InverseDisplacementField.IsNull()) {
			this->m_InverseInterpolator->SetInputImage(
					this->m_InverseDisplacementField);
		}
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetFixedParameters(
		const ParametersType & fixedParameters) {
	SizeValueType elems = NDimensions * (NDimensions + 3);

	if (fixedParameters.Size() != elems * 2) {
		itkExceptionMacro("The fixed parameters are not the right size.");
	}

	SizeType size;
	for (unsigned int d = 0; d < NDimensions; d++) {
		size[d] = static_cast<SizeValueType>(fixedParameters[d]);
	}

	PointType origin;
	for (unsigned int d = 0; d < NDimensions; d++) {
		origin[d] = fixedParameters[d + NDimensions];
	}

	SpacingType spacing;
	for (unsigned int d = 0; d < NDimensions; d++) {
		spacing[d] = fixedParameters[d + 2 * NDimensions];
	}

	DirectionType direction;
	for (unsigned int di = 0; di < NDimensions; di++) {
		for (unsigned int dj = 0; dj < NDimensions; dj++) {
			direction[di][dj] = fixedParameters[3 * NDimensions
					+ (di * NDimensions + dj)];
		}
	}

	PixelType zeroDisplacement;
	zeroDisplacement.Fill(0.0);

	typename DisplacementFieldType::Pointer displacementField =
			DisplacementFieldType::New();
	displacementField->SetSpacing(spacing);
	displacementField->SetOrigin(origin);
	displacementField->SetDirection(direction);
	displacementField->SetRegions(size);
	displacementField->Allocate();
	displacementField->FillBuffer(zeroDisplacement);

	this->SetDisplacementField(displacementField);

	if (!this->m_InverseDisplacementField.IsNull()) {
		typename DisplacementFieldType::Pointer inverseDisplacementField =
				DisplacementFieldType::New();
		inverseDisplacementField->SetSpacing(spacing);
		inverseDisplacementField->SetOrigin(origin);
		inverseDisplacementField->SetDirection(direction);
		inverseDisplacementField->SetRegions(size);
		inverseDisplacementField->Allocate();
		inverseDisplacementField->FillBuffer(zeroDisplacement);

		this->SetInverseDisplacementField(inverseDisplacementField);
	}

	// Set coefficients field
	for (unsigned int d = 0; d < NDimensions; d++) {
		size[d] = static_cast<SizeValueType>(fixedParameters[d + elems]);
	}

	for (unsigned int d = 0; d < NDimensions; d++) {
		origin[d] = fixedParameters[d + NDimensions + elems];
	}

	for (unsigned int d = 0; d < NDimensions; d++) {
		spacing[d] = fixedParameters[d + 2 * NDimensions + elems];
	}

	for (unsigned int di = 0; di < NDimensions; di++) {
		for (unsigned int dj = 0; dj < NDimensions; dj++) {
			direction[di][dj] = fixedParameters[3 * NDimensions
					+ (di * NDimensions + dj) + elems];
		}
	}

	typename DisplacementFieldType::Pointer coeffField =
			DisplacementFieldType::New();
	coeffField->SetSpacing(spacing);
	coeffField->SetOrigin(origin);
	coeffField->SetDirection(direction);
	coeffField->SetRegions(size);
	coeffField->Allocate();
	coeffField->FillBuffer(zeroDisplacement);

	this->SetCoefficientsField(coeffField);
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetFixedParametersFromDisplacementField() const {
	const typename DisplacementFieldType::RegionType & fieldRegion =
			this->m_DisplacementField->GetLargestPossibleRegion();

	// Set the field size parameters
	SizeType fieldSize = fieldRegion.GetSize();
	for (unsigned int i = 0; i < NDimensions; i++) {
		this->m_FixedParameters[i] =
				static_cast<ParametersValueType>(fieldSize[i]);
	}

	// Set the origin parameters
	PointType fieldOrigin = this->m_DisplacementField->GetOrigin();
	for (unsigned int i = 0; i < NDimensions; i++) {
		this->m_FixedParameters[NDimensions + i] = fieldOrigin[i];
	}

	// Set the spacing parameters
	SpacingType fieldSpacing = this->m_DisplacementField->GetSpacing();
	for (unsigned int i = 0; i < NDimensions; i++) {
		this->m_FixedParameters[2 * NDimensions + i] =
				static_cast<ParametersValueType>(fieldSpacing[i]);
	}

	// Set the direction parameters
	DirectionType fieldDirection = this->m_DisplacementField->GetDirection();
	for (unsigned int di = 0; di < NDimensions; di++) {
		for (unsigned int dj = 0; dj < NDimensions; dj++) {
			this->m_FixedParameters[3 * NDimensions + (di * NDimensions + dj)] =
					static_cast<ParametersValueType>(fieldDirection[di][dj]);
		}
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::SetFixedParametersFromCoefficientsField() const {
	SizeValueType elems = NDimensions * (NDimensions + 3);

	const typename DisplacementFieldType::RegionType & coeffRegion =
			this->m_CoefficientsField->GetLargestPossibleRegion();

	// Set the field size parameters
	SizeType coeffSize = coeffRegion.GetSize();
	for (unsigned int i = 0; i < NDimensions; i++) {
		this->m_FixedParameters[i + elems] =
				static_cast<ParametersValueType>(coeffSize[i]);
	}

	// Set the origin parameters
	PointType coeffOrigin = this->m_CoefficientsField->GetOrigin();
	for (unsigned int i = 0; i < NDimensions; i++) {
		this->m_FixedParameters[NDimensions + i + elems] = coeffOrigin[i];
	}

	// Set the spacing parameters
	SpacingType coeffSpacing = this->m_CoefficientsField->GetSpacing();
	for (unsigned int i = 0; i < NDimensions; i++) {
		this->m_FixedParameters[2 * NDimensions + i + elems] =
				static_cast<ParametersValueType>(coeffSpacing[i]);
	}

	// Set the direction parameters
	DirectionType coeffDirection = this->m_CoefficientsField->GetDirection();
	for (unsigned int di = 0; di < NDimensions; di++) {
		for (unsigned int dj = 0; dj < NDimensions; dj++) {
			this->m_FixedParameters[3 * NDimensions + (di * NDimensions + dj) + elems] =
					static_cast<ParametersValueType>(coeffDirection[di][dj]);
		}
	}
}

template<typename TScalar, unsigned int NDimensions>
void RBFFieldTransform<TScalar, NDimensions>::PrintSelf(
		std::ostream& os, Indent indent) const {
	Superclass::PrintSelf(os, indent);

	std::cout << indent << "Interpolator: " << std::endl;
	std::cout << indent << indent << this->m_Interpolator << std::endl;

	std::cout << indent << "InverseInterpolator: " << std::endl;
	std::cout << indent << indent << this->m_InverseInterpolator << std::endl;

	if (this->m_DisplacementField) {
		std::cout << indent << "Displacement Field: " << std::endl;
		std::cout << indent << indent << this->m_DisplacementField << std::endl;
	} else {
		std::cout << "Displacement field not set." << std::endl;
	}

	if (this->m_CoefficientsField) {
		std::cout << indent << "Coefficients Field: " << std::endl;
		std::cout << indent << indent << this->m_CoefficientsField << std::endl;
	} else {
		std::cout << "Coefficients field not set." << std::endl;
	}

	if (this->m_InverseDisplacementField) {
		std::cout << indent << "Inverse Displacement Field: " << std::endl;
		std::cout << indent << indent << this->m_InverseDisplacementField
				<< std::endl;
	} else {
		std::cout << "Inverse Displacement field not set." << std::endl;
	}
}

} // namespace itk

#endif
