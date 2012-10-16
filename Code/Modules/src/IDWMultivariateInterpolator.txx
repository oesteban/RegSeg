/* --------------------------------------------------------------------------------------
 * File:    IDWMultivariateInterpolator.txx
 * Date:    Mar 28, 2012
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2012, Oscar Esteban - oesteban@die.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM)
 and Signal Processing Lab 5, EPFL (LTS5-EPFL)
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names
 of its contributors may be used to endorse or promote products derived
 from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef IDWMULTIVARIATEINTERPOLATOR_TXX_
#define IDWMULTIVARIATEINTERPOLATOR_TXX_

#include "IDWMultivariateInterpolator.h"

namespace rstk {

template< class TInputMesh, class TOutputImage >
IDWMultivariateInterpolator< TInputMesh, TOutputImage >::
IDWMultivariateInterpolator() : m_pFactor(InputImageDimension) {
};


template< class TInputMesh, class TOutputImage >
typename IDWMultivariateInterpolator< TInputMesh, TOutputImage >::OutputType
IDWMultivariateInterpolator< TInputMesh, TOutputImage >::
Evaluate( const PointType & point ) const {
	// For every point in the mesh (c_i)
	PointType ci;
	PixelType vi, Tv;
	Tv.Fill(0.0);
	double wi, Tw = 0.0;

	for( size_t i = 0; i < this->m_NumberOfMeshPoints; i++ ) {
		// Get point
		ci = this->m_Mesh->GetPoint( i );
		this->m_Mesh->GetPointData( i, &vi );

		double dist = (ci - point).GetNorm();

		if( dist < vnl_math::eps ) {
			return vi;
		}

		// Compute weight i: wi(x) (inverse distance)
		wi = 1.0 / vcl_pow( dist, m_pFactor );

		// Accumulate weight
		Tw += wi;

		// Multiply by mesh's vector & Accumulate weighted vector
		Tv += (vi * wi);

	}
	// Divide accumulated weighted vector by accumulated total weight.
	return (Tv / Tw);



};

template< class TInputMesh, class TOutputImage >
void IDWMultivariateInterpolator< TInputMesh, TOutputImage >::
PrintSelf( std::ostream &os, Indent indent) const {
	  this->Superclass::PrintSelf(os, indent);
};

}

#endif /* IDWMULTIVARIATEINTERPOLATOR_TXX_ */