// --------------------------------------------------------------------------
// File:             MahalanobisLevelSets.hxx
// Date:             27/10/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWERegistration-Debug@Debug
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef MAHALANOBISLEVELSETS_HXX_
#define MAHALANOBISLEVELSETS_HXX_

#include "MahalanobisLevelSets.h"

namespace rstk {
template <class TTargetImage, class TDeformationField, class TContourDeformation>
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::MahalanobisLevelSets() {


}

template <class TTargetImage, class TDeformationField, class TContourDeformation>
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::ValueType
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::GetValue() const {
	// for all classes

		// if pixel is inside object

			// compute mahalanobis distance to class

			// add to energy

	// reset m_Vule and return
	return this->m_Value;
}


template <class TTargetImage, class TDeformationField, class TContourDeformation>
void
MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>
::GetLevelSetsMap( MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::ContourDeformationType & contourDeformation,
		           MahalanobisLevelSets<TTargetImage,TDeformationField,TContourDeformation>::DeformationFieldType & targetDeformation) const {
	// Compute mesh of normals
	NormalFilterPointer normals = NormalFilterType::New();
	normals->SetInput( contourDeformation );
	normals->Update();

	// for all node in contourDeformation

		// compute mahalanobis distance in position

	    // project to normal

	// Interpolate sparse velocity field to targetDeformation

}

}

#endif /* MAHALANOBISLEVELSETS_HXX_ */
