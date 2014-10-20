// --------------------------------------------------------------------------
// File:             DisplacementFieldFileWriter.h
// Date:             01/11/2012
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

#ifndef DISPLACEMENTFIELDFILEWRITER_H_
#define DISPLACEMENTFIELDFILEWRITER_H_

// Include headers
#include <itkObject.h>
#include <itkImageFileWriter.h>

// Namespace declaration

namespace rstk {
/** \class DisplacementFieldFileWriter
 *  \brief This class
 *
 *  Long description
 *
 *  \ingroup
 */

template< typename TDisplacementField >
class DisplacementFieldFileWriter: public itk::Object {
public:
	typedef DisplacementFieldFileWriter   Self;
	typedef itk::Object                   Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( DisplacementFieldFileWriter, itk::Object );
	itkNewMacro( Self );

	typedef TDisplacementField                              DisplacementFieldType;
	typedef typename DisplacementFieldType::ConstPointer    DisplacementFieldPointer;
	typedef typename DisplacementFieldType::PixelType       VectorType;
	typedef typename VectorType::ValueType                  ValueType;

	itkStaticConstMacro( Dimension, unsigned int, itkGetStaticConstMacro(VectorType::Dimension) );

	itkSetConstObjectMacro(Input,DisplacementFieldType);
	itkSetStringMacro(FileName);

	void Update() const {
		typedef itk::Image<ValueType,Dimension+1> FieldType;
		typename FieldType::Pointer out = FieldType::New();
		typename FieldType::SizeType outSize;
		typename FieldType::SpacingType outSpacing;
		for( size_t i = 0; i<Dimension; i++) {
			outSize[i] = m_Input->GetLargestPossibleRegion().GetSize()[i];
			outSpacing[i] = m_Input->GetSpacing()[i];
		}
		outSize[Dimension] = Dimension;
		outSpacing[Dimension] = 1.0;
		out->SetRegions( outSize );
		out->SetSpacing( outSpacing );
		out->Allocate();
		out->FillBuffer(0.0);

		ValueType* buffer = out->GetBufferPointer();
		const VectorType* vectBuffer = m_Input->GetBufferPointer();
		size_t nVect = m_Input->GetLargestPossibleRegion().GetNumberOfPixels();
		size_t nPix = out->GetLargestPossibleRegion().GetNumberOfPixels();
		size_t comp, idx2;
		for(size_t pix = 0; pix< nPix; pix++) {
			comp = pix / nVect;
			idx2 = pix % nVect;
			VectorType val = *(vectBuffer+idx2);
			*(buffer+pix) = val[comp];
		}

		typename itk::ImageFileWriter<FieldType>::Pointer w = itk::ImageFileWriter<FieldType>::New();
		w->SetInput( out );
		w->SetFileName( m_FileName );
		w->Update();
	}

protected:
	DisplacementFieldFileWriter(){}
	~DisplacementFieldFileWriter(){}

	void PrintSelf( std::ostream& os, itk::Indent indent) const {
		Superclass::Print(os, indent);
	}

private:
	DisplacementFieldFileWriter( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	std::string m_FileName;
	DisplacementFieldPointer m_Input;
}; // End of class DisplacementFieldFileWriter
} // End of namespace


#endif /* DISPLACEMENTFIELDFILEWRITER_H_ */
