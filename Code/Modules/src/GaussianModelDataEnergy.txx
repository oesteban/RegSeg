// --------------------------------------------------------------------------------------
// File:          GaussianModelDataEnergy.txx
// Date:          Apr 3, 2012
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       0.1
// License:       BSD
// Short Summary: 
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Biomedical Image Technology, UPM (BIT-UPM)
// and Signal Processing Lab 5, EPFL (LTS5-EPFL)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//   * Neither the names of the copyright holders, nor the names of
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
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

#ifndef GAUSSIANMODELDATAENERGY_TXX_
#define GAUSSIANMODELDATAENERGY_TXX_


namespace rstk {

template< class TFixedImage, class TMovingSurface, class TEnergyValue >
GaussianModelDataEnergy< TFixedImage, TMovingSurface, TEnergyValue >::
GaussianModelDataEnergy() {
	m_Interpolator = DefaultInterpolatorType::New();
}

template< class TFixedImage, class TMovingSurface, class TEnergyValue >
typename GaussianModelDataEnergy< TFixedImage, TMovingSurface, TEnergyValue >::SparseVectorFieldType*
GaussianModelDataEnergy< TFixedImage, TMovingSurface, TEnergyValue >::
ComputeNormals() {

}


template< class TFixedImage, class TMovingSurface, class TEnergyValue >
typename GaussianModelDataEnergy< TFixedImage, TMovingSurface, TEnergyValue >::EnergyValueType
GaussianModelDataEnergy< TFixedImage, TMovingSurface, TEnergyValue >::
GetGradient() {
	EnergyValueType gradient, Ed = 0.0;

	MovingSurfaceConstPointer surf = m_DisplacementField->GetTransformedSurface();
	NormalsCalculatorPointer n = NormalsCalculatorFilter::New();
	n->SetInput( surf );
	n->Update();

	NormalsMeshConstPointer normalSurf = n->GetOutput();
	FeatureVectorType feature;

	size_t nPoints = surf->GetNumberOfPoints();
	PointType p;
	DisplacementVectorType normal;
	DisplacementVectorType wi;
	vnl_vector< double > tempVector(FeatureVectorSize);

	for ( size_t i = 0; i < nPoints; i++) {           // For every point in exterior mesh
		// Get Point position
		surf->GetPoint( i, &p );

		// Get interpolated feature vector at pixel position p
		feature = m_FixedImage->GetPixel( );

		// Get normal of surface in that point
		normalSurf->GetPointData(i, &normal);

		// Get k-th weight (RBF interpolator)
		wi = m_Interpolator->Evaluate( p );


		// Compute energy change with respect to the model
		for ( size_t i = 0; i < FeatureVectorSize; ++i ) {
			tempVector[i] = feature[i] - m_ModelMean[i];
		}
		Ed = dot_product( tempVector, m_ModelInvCovariance.GetVnlMatrix() * tempVector )
				* (normal*wi);
		gradient += Ed;
	}


	return (-gradient);
}

}


#endif /* GAUSSIANMODELDATAENERGY_TXX_ */
