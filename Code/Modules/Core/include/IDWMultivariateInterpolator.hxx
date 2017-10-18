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

#ifndef IDWMULTIVARIATEINTERPOLATOR_HXX_
#define IDWMULTIVARIATEINTERPOLATOR_HXX_


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


#endif /* IDWMULTIVARIATEINTERPOLATOR_HXX_ */
