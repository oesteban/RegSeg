// --------------------------------------------------------------------------
// File:             ContourDisplacementField.h
// Date:             29/10/2012
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

#ifndef CONTOURDISPLACEMENTFIELD_H_
#define CONTOURDISPLACEMENTFIELD_H_

// Include headers
#include <itkObject.h>
#include <itkVector.h>
#include <itkMesh.h>
#include <itkQuadEdgeMesh.h>

// Namespace declaration

namespace rstk {
/** \class ContourDisplacementField
 *  \brief This class implements sparse displacement field defined along a reference mesh (discrete 2D-manifold)
 *
 *  ContourDisplacementField implements a QuadEdgeMesh in order to having available all the QE filters (specifically,
 *  the Normals filter). It is built upon a FreeSurfer's mesh set via SetReferenceMesh(). ContourDisplacementField
 *  process 3D deformation surfaces and fields by default.
 *
 *  \ingroup RSTKCore
 */

template< typename TCoordPrecision = double, unsigned int VDimension = 3u >
class ContourDisplacementField: public itk::QuadEdgeMesh< itk::Vector< TCoordPrecision, VDimension >, VDimension > {
public:
	typedef ContourDisplacementField      Self;
	typedef itk::QuadEdgeMesh< itk::Vector< TCoordPrecision, VDimension >, VDimension >
	                                      Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( ContourDisplacementField, itk::Object );
	itkNewMacro( Self );

	typedef typename Superclass::PixelType                      PixelType;
	typedef itk::Mesh<TCoordPrecision, VDimension >             InternalMeshType;
	typedef typename InternalMeshType::Pointer                  InternalMeshPointer;

	itkSetObjectMacro( ReferenceMesh, InternalMeshType );
	itkGetConstObjectMacro( ReferenceMesh, InternalMeshType );

protected:
	ContourDisplacementField();
	virtual ~ContourDisplacementField() {}

	void PrintSelf( std::ostream& os, itk::Indent indent) const;

private:
	ContourDisplacementField( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	InternalMeshPointer m_ReferenceMesh;

}; // End of class ContourDisplacementField
} // End of namespace


#endif /* CONTOURDISPLACEMENTFIELD_H_ */
