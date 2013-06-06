// --------------------------------------------------------------------------
// File:             MeshToImageFilter.h
// Date:             30/10/2012
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
// This file is part of ACWEReg
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

#ifndef MESHTOIMAGEFILTER_H_
#define MESHTOIMAGEFILTER_H_


#include "itkImageSource.h"


using namespace itk;

namespace rstk {

/** \class MeshToImageFilter
 * \brief
 *
 * MeshToImageFilter is the base class for all process objects that output
 * image data and require a Mesh data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ImageFilters
 * \ingroup ITKMesh
 */
template< class TInputMesh, class TOutputImage >
class MeshToImageFilter:public ImageSource< TOutputImage >
{
public:
	typedef MeshToImageFilter                              Self;
	typedef ImageSource< TOutputImage >                    Superclass;
	typedef SmartPointer< Self >                           Pointer;
	typedef SmartPointer< const Self >                     ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( MeshToImageFilter, ImageSource );

	typedef          TInputMesh                            InputMeshType;
	typedef typename InputMeshType::Pointer                InputMeshPointer;
	typedef typename InputMeshType::ConstPointer           InputMeshConstPointer;
	typedef typename InputMeshType::PointType              MeshPointType;
	typedef          TOutputImage                          OutputImageType;
	typedef typename OutputImageType::Pointer              OutputImagePointer;
	typedef typename OutputImageType::ConstPointer         OutputImageConstPointer;
	typedef typename OutputImageType::SizeType             OutputImageSizeType;
	typedef typename OutputImageType::PixelType            PointType;

	typedef typename std::vector< InputMeshConstPointer >  InputMeshList;

	itkStaticConstMacro( Dimension, unsigned int, TOutputImage::ImageDimension );


	virtual void SetInput( const InputMeshType *input ) {
		this->SetInput( 0, input );
	}
	virtual void SetInput( size_t idx, const InputMeshType *input) {
		if ( idx >= m_Mesh.size() ) m_Mesh.resize( idx+1 );

		this->m_Mesh[idx] = input;
	}

	const InputMeshType*   GetInput(void) { return this->GetInput(0); }
	const InputMeshType*   GetInput(unsigned int idx) { return this->m_Mesh[idx]; }

	size_t GetNumberOfIndexedInputs() { return m_Mesh.size(); }

protected:
	MeshToImageFilter(){};
	~MeshToImageFilter(){};

	//virtual void GenerateData();
	//virtual void PrintSelf( std::ostream &os, Indent indent) const;

	InputMeshList m_Mesh;

private:
	MeshToImageFilter( const Self& ); // purposely not implemented
	void operator=( const Self& );    // purposely not implemented


};

}


#endif /* MESHTOIMAGEFILTER_H_ */
