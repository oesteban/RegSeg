// --------------------------------------------------------------------------
// File:             ComponentsFileWriter.h
// Date:             22/07/2013
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban)
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

#ifndef ComponentsFileWriter_H_
#define ComponentsFileWriter_H_

// Include headers
#include <itkObject.h>
#include <itkImageFileWriter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>

#include <iostream>
// Namespace declaration

namespace rstk {
/** \class ComponentsFileWriter
 *  \brief This class
 *
 *  Long description
 *
 *  \ingroup
 */

template< typename TVectorImageType >
class ComponentsFileWriter: public itk::Object {
public:
	typedef ComponentsFileWriter   Self;
	typedef itk::Object                   Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( ComponentsFileWriter, itk::Object );
	itkNewMacro( Self );

	typedef TVectorImageType                              VectorImageType;
	typedef typename VectorImageType::Pointer             VectorImagePointer;
	typedef typename VectorImageType::ConstPointer        VectorImageConstPointer;
	typedef typename VectorImageType::PixelType           VectorType;
	typedef typename VectorType::ValueType                ValueType;

	itkStaticConstMacro( Dimension, unsigned int, itkGetStaticConstMacro(VectorImageType::ImageDimension) );

	typedef typename itk::Image< ValueType, Dimension >     ComponentImageType;
	typedef itk::ImageFileWriter<ComponentImageType >       ComponentWriterType;
	typedef typename ComponentWriterType::Pointer           ComponentWriterPointer;

	typedef itk::VectorIndexSelectionCastImageFilter< VectorImageType, ComponentImageType >
	                                                        SelectComponentFilter;
	typedef typename SelectComponentFilter::Pointer         SelectComponentPointer;

	itkSetConstObjectMacro(Input, VectorImageType);
	itkSetStringMacro(FileName);

	void Update() const {
		size_t ncomps = this->m_Input->GetNumberOfComponentsPerPixel();

		SelectComponentPointer adaptor = SelectComponentFilter::New();
		adaptor->SetInput(this->m_Input);

		for( size_t comp=0; comp<ncomps; comp++) {
			adaptor->SetIndex(comp);
			adaptor->Update();

			std::stringstream ss;
			ss << m_FileName << "_cmp" << comp << ".nii.gz";
			ComponentWriterPointer w = ComponentWriterType::New();
			w->SetInput( adaptor->GetOutput() );
			w->SetFileName( ss.str().c_str() );
			w->Update();
		}
	}

protected:
	ComponentsFileWriter(){}
	~ComponentsFileWriter(){}

	void PrintSelf( std::ostream& os, itk::Indent indent) const {
		Superclass::Print(os, indent);
	}

private:
	ComponentsFileWriter( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	std::string m_FileName;
	VectorImageConstPointer m_Input;
}; // End of class ComponentsFileWriter
} // End of namespace

#endif /* ComponentsFileWriter_H_ */
