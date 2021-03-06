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

#ifndef __MahalanobisDistanceMembershipFunction_hxx
#define __MahalanobisDistanceMembershipFunction_hxx

#include "MahalanobisDistanceMembershipFunction.h"

#include <itkLightObject.h>
#include <math.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>

namespace rstk {
template<typename TVector>
MahalanobisDistanceMembershipFunction<TVector>::MahalanobisDistanceMembershipFunction() {
	itk::NumericTraits<MeanVectorType>::SetLength(m_Mean,
			this->GetMeasurementVectorSize());
	m_Mean.Fill(0.0f);

	itk::NumericTraits<MeanVectorType>::SetLength(m_RangeMin,
			this->GetMeasurementVectorSize());
	m_RangeMin.Fill(itk::NumericTraits<MeasurementValueType>::min());

	itk::NumericTraits<MeanVectorType>::SetLength(m_RangeMax,
			this->GetMeasurementVectorSize());
	m_RangeMax.Fill(itk::NumericTraits<MeasurementValueType>::max());

	m_Covariance.SetSize(this->GetMeasurementVectorSize(),
			this->GetMeasurementVectorSize());
	m_Covariance.SetIdentity();
	m_InverseCovariance = m_Covariance;
	m_CovarianceNonsingular = true;

	m_MaximumValue = itk::NumericTraits<double>::max();
}

template<typename TVector>
void MahalanobisDistanceMembershipFunction<TVector>::SetMean(
		const MeanVectorType & mean) {
	if (this->GetMeasurementVectorSize()) {
		itk::Statistics::MeasurementVectorTraits::Assert(mean,
				this->GetMeasurementVectorSize(),
				"GaussianMembershipFunction::SetMean(): Size of mean vector specified does not match the size of a measurement vector.");
	} else {
		// not already set, cache the size
		this->SetMeasurementVectorSize(mean.Size());
	}

	if (m_Mean != mean) {
		m_Mean = mean;
		this->Modified();
	}
}

template<typename TVector>
void MahalanobisDistanceMembershipFunction<TVector>::SetRange(
		const MeanVectorType & lower, const MeanVectorType & upper) {
	if (this->GetMeasurementVectorSize()) {
		itk::Statistics::MeasurementVectorTraits::Assert(lower,
				this->GetMeasurementVectorSize(),
				"MembershipFunction::SetRange(): Size of lower bound vector specified does not match the size of a measurement vector.");
		itk::Statistics::MeasurementVectorTraits::Assert(upper,
				this->GetMeasurementVectorSize(),
				"MembershipFunction::SetRange(): Size of upper bound vector specified does not match the size of a measurement vector.");
	} else {
		// not already set, cache the size
		this->SetMeasurementVectorSize(lower.Size());
	}

	if (m_RangeMax != upper) {
		m_RangeMax = upper;
		this->Modified();
	}
	if (m_RangeMin != lower) {
		m_RangeMin = lower;
		this->Modified();
	}
}

template<typename TVector>
void MahalanobisDistanceMembershipFunction<TVector>::SetCovariance(
		const CovarianceMatrixType & cov) {
	// Sanity check
	if (cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols()) {
		itkExceptionMacro(<< "Covariance matrix must be square");
	}
	if (this->GetMeasurementVectorSize()) {
		if (cov.GetVnlMatrix().rows() != this->GetMeasurementVectorSize()) {
			itkExceptionMacro(
					<< "Length of measurement vectors must be" << " the same as the size of the covariance.");
		}
	} else {
		// not already set, cache the size
		this->SetMeasurementVectorSize(cov.GetVnlMatrix().rows());
	}

	if (m_Covariance == cov) {
		// no need to copy the matrix, compute the inverse, or the normalization
		return;
	}

	m_Covariance = cov;

	// the inverse of the covariance matrix is first computed by SVD
	vnl_matrix_inverse<double> inv_cov(m_Covariance.GetVnlMatrix());

	// the determinant is then costless this way
	double det = inv_cov.determinant_magnitude();

	if (det < 0.) {
		itkExceptionMacro(<< "det( m_Covariance ) < 0");
	}

	m_OffsetTerm = log(2.0 * vnl_math::pi) * m_Covariance.Rows() + log(det);

	// 1e-6 is an arbitrary value!!!
	const double singularThreshold = 1.0e-6;
	m_CovarianceNonsingular = (det > singularThreshold);

	if (m_CovarianceNonsingular) {
		// allocate the memory for m_InverseCovariance matrix
		m_InverseCovariance.GetVnlMatrix() = inv_cov.inverse();
	} else {
		// define the inverse to be diagonal with large values along the
		// diagonal. value chosen so (X-M)'inv(C)*(X-M) will usually stay
		// below NumericTraits<double>::max()
		const double aLargeDouble = std::pow(itk::NumericTraits<double>::max(),
				1.0 / 3.0) / (double) this->GetMeasurementVectorSize();
		m_InverseCovariance.SetSize(this->GetMeasurementVectorSize(),
				this->GetMeasurementVectorSize());
		m_InverseCovariance.SetIdentity();
		m_InverseCovariance *= aLargeDouble;
	}

	this->Modified();
}

template<typename TVector>
void MahalanobisDistanceMembershipFunction<TVector>::Initialize() {
	double uppertmp = Evaluate(m_RangeMax);
	double lowertmp = Evaluate(m_RangeMin);

	m_MaximumValue = (uppertmp > lowertmp)?uppertmp:lowertmp;
}

template<typename TVector>
double MahalanobisDistanceMembershipFunction<TVector>::Evaluate(
		const MeasurementVectorType & measurement) const {
	const MeasurementVectorSizeType measurementVectorSize =
			this->GetMeasurementVectorSize();

	// Our inverse covariance is always well formed. When the covariance
	// is singular, we use a diagonal inverse covariance with a large diagnonal

	// temp = ( y - mean )^t * InverseCovariance * ( y - mean )
	//
	// This is manually done to remove dynamic memory allocation:
	// double temp = dot_product( tempVector,  m_InverseCovariance.GetVnlMatrix() * tempVector );
	//
	double temp = 0.0;
	for (unsigned int r = 0; r < measurementVectorSize; ++r) {
		double rowdot = 0.0;
		for (unsigned int c = 0; c < measurementVectorSize; ++c) {
			rowdot += m_InverseCovariance(r, c) * (measurement[c] - m_Mean[c]);
		}
		temp += rowdot * (measurement[r] - m_Mean[r]);

		if (temp > m_MaximumValue)
			return m_MaximumValue;
	}

	return temp;
}

template<typename TVector>
void MahalanobisDistanceMembershipFunction<TVector>::PrintSelf(
		std::ostream & os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "Mean: " << m_Mean << std::endl;
	os << indent << "Covariance: " << std::endl;
	os << m_Covariance.GetVnlMatrix();
	os << indent << "InverseCovariance: " << std::endl;
	os << indent << m_InverseCovariance.GetVnlMatrix();
	os << indent << "Covariance nonsingular: "
			<< (m_CovarianceNonsingular ? "true" : "false") << std::endl;
}

template<typename TVector>
typename itk::LightObject::Pointer MahalanobisDistanceMembershipFunction<TVector>::InternalClone() const {
	itk::LightObject::Pointer loPtr = Superclass::InternalClone();
	typename Self::Pointer membershipFunction =
			dynamic_cast<Self *>(loPtr.GetPointer());
	if (membershipFunction.IsNull()) {
		itkExceptionMacro(
				<< "downcast to type " << this->GetNameOfClass() << " failed.");
	}

	membershipFunction->SetMeasurementVectorSize(
			this->GetMeasurementVectorSize());
	membershipFunction->SetMean(this->GetMean());
	membershipFunction->SetCovariance(this->GetCovariance());
	membershipFunction->SetRange(this->GetRangeMin(), this->GetRangeMax());

	return loPtr;
}

template<typename TVector>
void MahalanobisDistanceMembershipFunction<TVector>::UpdateVectorSize() {
	MeasurementVectorType m;

	if (itk::Statistics::MeasurementVectorTraits::IsResizable(m)) {
		if (itk::NumericTraits<MeanVectorType>::GetLength(m_Mean)
				!= this->GetMeasurementVectorSize()) {
			itk::NumericTraits<MeanVectorType>::SetLength(m_Mean,
					this->GetMeasurementVectorSize());
			m_Mean.Fill(0.0f);
			this->Modified();
		}

		if (itk::NumericTraits<MeanVectorType>::GetLength(m_RangeMin)
				!= this->GetMeasurementVectorSize()) {
			itk::NumericTraits<MeanVectorType>::SetLength(m_RangeMin,
					this->GetMeasurementVectorSize());
			m_RangeMin.Fill(itk::NumericTraits<MeasurementValueType>::min());
			this->Modified();
		}

		if (itk::NumericTraits<MeanVectorType>::GetLength(m_RangeMax)
				!= this->GetMeasurementVectorSize()) {
			itk::NumericTraits<MeanVectorType>::SetLength(m_RangeMax,
					this->GetMeasurementVectorSize());
			m_RangeMax.Fill(itk::NumericTraits<MeasurementValueType>::max());
			this->Modified();
		}

		if (m_Covariance.Rows() != this->GetMeasurementVectorSize()
				|| m_Covariance.Cols() != this->GetMeasurementVectorSize()) {
			m_Covariance.SetSize(this->GetMeasurementVectorSize(),
					this->GetMeasurementVectorSize());
			m_Covariance.SetIdentity();
			m_InverseCovariance = m_Covariance;
			m_CovarianceNonsingular = true;
			this->Modified();
		}
	}
}
} // end of namespace itk

#endif
