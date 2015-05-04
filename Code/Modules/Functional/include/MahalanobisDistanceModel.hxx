// --------------------------------------------------------------------------------------
// File:          MahalanobisDistanceModel.hxx
// Date:          Dec 23, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
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
// along with ACWEReg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef _MAHALANOBISDISTANCEMODEL_HXX_
#define _MAHALANOBISDISTANCEMODEL_HXX_

#include "MahalanobisDistanceModel.h"


#include <math.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_ldl_cholesky.h>

#include <string>
#include <fstream>
#include <streambuf>
#include <jsoncpp/json/json.h>

#include <itkNumericTraitsVariableLengthVectorPixel.h>

namespace rstk {
template< typename TInputVectorImage, typename TPriorsPrecisionType >
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::MahalanobisDistanceModel(): Superclass() {}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::PrintSelf(std::ostream & os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::InitializeMemberships() {
	this->m_NumberOfRegions = this->GetPriorsMap()->GetNumberOfComponentsPerPixel();

	size_t nregions = this->m_NumberOfRegions - this->m_NumberOfSpecialRegions;
	this->m_Memberships.resize(this->m_NumberOfRegions);
	this->m_Means.resize(nregions);
	this->m_RangeLower.resize(nregions);
	this->m_RangeUpper.resize(nregions);
	this->m_Covariances.resize(nregions);
	this->m_RegionOffsetContainer.SetSize(nregions);
	this->m_RegionOffsetContainer.Fill(0.0);
}


template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::PostGenerateData() {
	// Append special memberships
	if (this->m_NumberOfSpecialRegions >= 2) {
		MeasureType lastMax = this->m_RegionOffsetContainer[this->m_NumberOfRegions - this->m_NumberOfSpecialRegions - 1];
		UniformFunctionPointer bg_mf = UniformFunctionType::New();
		bg_mf->SetValue(lastMax * 0.25);
		this->m_Memberships[this->m_NumberOfRegions - 2] = bg_mf;
	}

	if (this->m_NumberOfSpecialRegions >= 1) {
		UniformFunctionPointer offmask_mf = UniformFunctionType::New();
		offmask_mf->SetValue(this->m_MaxEnergy);
		this->m_Memberships[this->m_NumberOfRegions - 1] = offmask_mf;
	}

	size_t ncomps = itk::NumericTraits<MeasurementVectorType>::GetLength(this->m_Means[0]);
	itk::NumericTraits<MeasurementVectorType>::SetLength(m_InvalidValue, ncomps);
	for( size_t i = 0; i < ncomps; i++ ) {
		m_InvalidValue = 0.0;
	}
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::GenerateData() {
	this->InitializeMemberships();
	this->EstimateRobust();
	this->PostGenerateData();
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::Estimate() {
	ReferenceSamplePointer sample = ReferenceSampleType::New();
	sample->SetImage( this->GetInput() );
	size_t npix = sample->Size();
	size_t nregions = this->m_NumberOfRegions - this->m_NumberOfSpecialRegions;

	const PriorsPrecisionType* priors = this->GetPriorsMap()->GetBufferPointer();
	size_t offset = this->GetPriorsMap()->GetNumberOfComponentsPerPixel();

	std::vector<WeightArrayType> weights;
	for( size_t roi = 0; roi < nregions; roi++ ) {
		WeightArrayType warr;
		warr.SetSize(npix);
		warr.Fill(0.0);
		weights.push_back(warr);
	}

	PriorsPrecisionType w;
	for( size_t i = 0; i < npix; i++ ) {
		for( size_t roi = 0; roi < nregions; roi++ ) {
			w = *(priors + offset * i + roi);
			weights[roi][i] = w;
		}
	}

	for( size_t roi = 0; roi < nregions; roi++ ) {
		CovarianceFilterPointer covFilter = CovarianceFilter::New();
		covFilter->SetInput( sample );
		covFilter->SetWeights( weights[roi] );
		covFilter->Update();

		InternalFunctionPointer mf = InternalFunctionType::New();
		mf->SetMean( covFilter->GetMean() );

		CovarianceMatrixType cov = covFilter->GetCovarianceMatrix();
		mf->SetCovariance( cov );

		this->m_Memberships[roi] = mf;
		this->m_RegionOffsetContainer[roi] = log(2.0 * vnl_math::pi) * cov.Rows() + this->ComputeCovarianceDeterminant(cov);
		this->m_Means[roi] = covFilter->GetMean();
		this->m_Covariances[roi] = cov;
	}
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::EstimateRobust() {
	size_t npix = this->GetInput()->GetLargestPossibleRegion().GetNumberOfPixels();
	size_t ncomps = this->GetInput()->GetNumberOfComponentsPerPixel();

	size_t offset = this->GetPriorsMap()->GetNumberOfComponentsPerPixel();
	size_t nregions = this->m_NumberOfRegions - this->m_NumberOfSpecialRegions;

	const PriorsPrecisionType* priors = this->GetPriorsMap()->GetBufferPointer();
	const PixelValueType * input = this->GetInput()->GetBufferPointer();


	std::vector<WeightArrayType> weights;
	for( size_t roi = 0; roi < nregions; roi++ ) {
		WeightArrayType warr;
		warr.SetSize(npix);
		warr.Fill(0.0);
		weights.push_back(warr);
	}

	PriorsPrecisionType w;
	for( size_t i = 0; i < npix; i++ ) {
		if( *(input + i * ncomps) == 0 && *(input + i * ncomps + 1) == 0) {
			continue;
		}

		for( size_t roi = 0; roi < nregions; roi++ ) {
			w = *(priors + offset * i + roi);
			weights[roi][i] = w;
		}
	}

	ReferenceSamplePointer sample = ReferenceSampleType::New();
	sample->SetImage( this->GetInput() );
	for( size_t roi = 0; roi < nregions; roi++ ) {
		InternalFunctionPointer mf = InternalFunctionType::New();

		CovarianceFilterPointer covFilter = CovarianceFilter::New();
		covFilter->SetInput( sample );
		covFilter->SetWeights( weights[roi] );
		covFilter->Update();

		MeasurementVectorType mean = covFilter->GetMean();
		mf->SetMean(mean);

		CovarianceMatrixType cov = covFilter->GetCovarianceMatrix();
		mf->SetCovariance( cov );

		MeasurementVectorType range[2];

		this->m_RangeLower[roi] = covFilter->GetRangeMin();
		this->m_RangeUpper[roi] = covFilter->GetRangeMax();
		mf->SetRange(this->m_RangeLower[roi], this->m_RangeUpper[roi]);

		mf->Initialize();

		this->m_Memberships[roi] = mf;
		this->m_RegionOffsetContainer[roi] = mf->GetOffsetTerm();

		double maxv = mf->GetMaximumValue() * 1.0e3;
		if( maxv > this->m_MaxEnergy ) {
			this->m_MaxEnergy = maxv;
		}
		this->m_Means[roi] = mean;
		this->m_Covariances[roi] = cov;
	}
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
typename MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >::MeasureType
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::ComputeCovarianceDeterminant(CovarianceMatrixType& cov) const{
	typedef typename CovarianceMatrixType::ValueType ValueType;
	size_t ncomps = cov.Rows();
	CovarianceMatrixType invcov;
	MeasureType det = 1.0;

	if( ncomps > 1 ) {
		// Compute diagonal and check that eigenvectors >= 0.0
		typedef typename vnl_diag_matrix<ValueType>::iterator DiagonalIterator;
		typedef vnl_symmetric_eigensystem<ValueType> Eigensystem;
		vnl_matrix< ValueType > vnlCov = cov.GetVnlMatrix();
		Eigensystem* e = new Eigensystem( vnlCov );

		bool modified = false;
		DiagonalIterator itD = e->D.begin();
		while ( itD!= e->D.end() ) {
			if (*itD < 0) {
				*itD = 0.;
				modified = true;
			}
			itD++;
		}

		if (modified)
			cov = e->recompose();

		delete e;

		// the inverse of the covariance matrix is first computed by SVD
		vnl_matrix_inverse< ValueType > inv_cov( cov.GetVnlMatrix() );

		// the determinant is then costless this way
		det = inv_cov.determinant_magnitude();

		if( det < 0.) {
			itkExceptionMacro( << "| sigma | < 0" );
		}

		// FIXME Singurality Threshold for Covariance matrix: 1e-6 is an arbitrary value!!!
		const ValueType singularThreshold = 1.0e-10;
		if( det > singularThreshold ) {
			// allocate the memory for inverse covariance matrix
			invcov = inv_cov.inverse();
		} else {
			// TODO Perform cholesky diagonalization and select the semi-positive aproximation
			vnl_matrix< double > diag_cov( ncomps, ncomps );
			for ( size_t i = 0; i<ncomps; i++)
				for ( size_t j = 0; j<ncomps; j++)
					diag_cov[i][j] = vnlCov[i][j];
			vnl_ldl_cholesky* chol = new vnl_ldl_cholesky( diag_cov );
			vnl_vector< double > D( chol->diagonal() );
			det = dot_product( D, D );
			vnl_matrix_inverse< double > R (chol->upper_triangle());

			for ( size_t i = 0; i<ncomps; i++)
				for ( size_t j = 0; j<ncomps; j++)
					invcov(i,j) = R.inverse()[i][j];
		}
	} else {
		det = fabs( cov(0,0) );
	}
	return log( det );
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
std::string
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::PrintFormattedDescriptors() {
	Json::Value root = Json::Value( Json::objectValue );

	size_t nrois = this->m_NumberOfRegions - this->m_NumberOfSpecialRegions;
	root["descriptors"]["number"] = static_cast<Json::UInt>(nrois);
	root["descriptors"]["special"] = static_cast<Json::UInt>(this->m_NumberOfSpecialRegions);
	root["descriptors"]["values"] = Json::Value( Json::arrayValue );
	root["descriptors"]["max_energy"] = Json::Value( this->m_MaxEnergy );

	for ( size_t i = 0; i<nrois; i++ ){
		Json::Value vnode = Json::Value( Json::objectValue );
		vnode["id"] = static_cast<Json::UInt>( i );
		//vnode["determinant"] = Json::Value(this->ComputeCovarianceDeterminant(this->m_Covariances[i]));
		vnode["offset"] = Json::Value( this->m_RegionOffsetContainer[i] );
		vnode["mu"] = Json::Value( Json::arrayValue );
		vnode["cov"] = Json::Value( Json::arrayValue );
		vnode["range"]["lower"] = Json::Value( Json::arrayValue );
		vnode["range"]["upper"] = Json::Value( Json::arrayValue );
		size_t ncomps = itk::NumericTraits<MeasurementVectorType>::GetLength(this->m_Means[i]);

		for (size_t k = 0; k < ncomps; k++) {
			vnode["mu"].append(Json::Value(this->m_Means[i][k]));
			vnode["range"]["lower"].append(Json::Value(this->m_RangeLower[i][k]));
			vnode["range"]["upper"].append(Json::Value(this->m_RangeUpper[i][k]));

			for (size_t j = 0; j < ncomps; j++) {
				vnode["cov"].append(Json::Value(this->m_Covariances[i][j][k]));
			}

		}

		root["descriptors"]["values"].append(vnode);
	}

	return root.toStyledString();
}

template< typename TInputVectorImage, typename TPriorsPrecisionType >
void
MahalanobisDistanceModel< TInputVectorImage, TPriorsPrecisionType >
::ReadDescriptorsFromFile(std::string filename) {
	Json::Reader reader;
	std::ifstream t(filename.c_str());
	std::string str((std::istreambuf_iterator<char>(t)),
			         std::istreambuf_iterator<char>());

	Json::Value root;
	bool parsed = reader.parse(str, root, false);

	if (!parsed) {
		itkExceptionMacro(<< "Failed to read JSON file: " << reader.getFormattedErrorMessages());
	}

	const size_t nregions = root["descriptors"]["number"].asUInt();
	this->m_NumberOfSpecialRegions = root["descriptors"]["special"].asUInt();
	this->m_NumberOfRegions = nregions + this->m_NumberOfSpecialRegions;

	const Json::Value values = root["descriptors"]["values"];

	if (nregions != values.size()) {
		itkExceptionMacro(<< "Failed to read JSON file: number of regions does not match.");
	}

	this->m_Memberships.resize(this->m_NumberOfRegions);
	this->m_Means.resize(nregions);
	this->m_RangeLower.resize(nregions);
	this->m_RangeUpper.resize(nregions);
	this->m_Covariances.resize(nregions);
	this->m_RegionOffsetContainer.SetSize(nregions);
	this->m_RegionOffsetContainer.Fill(0.0);

	for(unsigned int roi = 0; roi < nregions; roi++ ) {
    	const Json::Value regnode = values[roi];
    	size_t ncomps = regnode["mu"].size();

    	MeasurementVectorType mean;
    	MeasurementVectorType lower;
    	MeasurementVectorType upper;
    	CovarianceMatrixType cov(ncomps, ncomps);
    	mean.SetSize(ncomps);
    	lower.SetSize(ncomps);
    	upper.SetSize(ncomps);
    	for(unsigned int i = 0; i < ncomps; i++) {
    		mean[i] = regnode["mu"][i].asFloat();

    		lower[i] = regnode["range"]["lower"][i].asFloat();
    		upper[i] = regnode["range"]["upper"][i].asFloat();

    		for(size_t j = 0; j < ncomps; j++) {
    			unsigned int cv = i * ncomps + j;
    			cov(j, i) = regnode["cov"][cv].asFloat();
    		}
    	}

    	this->m_RangeLower[roi] = lower;
    	this->m_RangeUpper[roi] = upper;
		this->m_Means[roi] = mean;
		this->m_Covariances[roi] = cov;
		this->m_RegionOffsetContainer[roi] = regnode["offset"].asFloat();

		InternalFunctionPointer mf = InternalFunctionType::New();
		mf->SetMean(mean);
		mf->SetCovariance( cov );
		mf->SetRange( lower, upper );
		mf->Initialize();

		this->m_Memberships[roi] = mf;
	}
	this->m_MaxEnergy = root["descriptors"]["max_energy"].asFloat();

	this->PostGenerateData();
}

}


#endif /* _MAHALANOBISDISTANCEMODEL_HXX_ */
