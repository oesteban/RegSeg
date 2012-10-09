/* --------------------------------------------------------------------------------------
 * File:    MeshToImageFilter.h
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
	typedef          TOutputImage                          OutputImageType;
	typedef typename OutputImageType::Pointer              OutputImagePointer;
	typedef typename OutputImageType::ConstPointer         OutputImageConstPointer;
	typedef typename OutputImageType::SizeType             OutputImageSizeType;
	typedef typename OutputImageType::PixelType            ValueType;

	itkStaticConstMacro( OutputImageDimension, unsigned int, TOutputImage::ImageDimension );


	void SetInput( const InputMeshType *input ) { m_Mesh = input; }
	const InputMeshType* GetInput(void) { return m_Mesh; }

protected:
	MeshToImageFilter(){};
	~MeshToImageFilter(){};

	//virtual void GenerateData();
	//virtual void PrintSelf( std::ostream &os, Indent indent) const;

	InputMeshConstPointer m_Mesh;

private:
	MeshToImageFilter( const Self& ); // purposely not implemented
	void operator=( const Self& );    // purposely not implemented


};

}

#endif /* MESHTOIMAGEFILTER_H_ */
