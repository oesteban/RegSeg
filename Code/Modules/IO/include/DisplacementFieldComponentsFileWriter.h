// --------------------------------------------------------------------------
// File:             DisplacementFieldComponentsFileWriter.h
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

#ifndef DisplacementFieldComponentsFileWriter_H_
#define DisplacementFieldComponentsFileWriter_H_

// Include headers
#include <itkObject.h>
#include <itkImageFileWriter.h>

#include <iostream>
// Namespace declaration

namespace rstk {
/** \class DisplacementFieldComponentsFileWriter
 *  \brief This class
 *
 *  Long description
 *
 *  \ingroup
 */

template< typename TDisplacementField >
class DisplacementFieldComponentsFileWriter: public itk::Object {
public:
	typedef DisplacementFieldComponentsFileWriter   Self;
	typedef itk::Object                   Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( DisplacementFieldComponentsFileWriter, itk::Object );
	itkNewMacro( Self );

	typedef TDisplacementField                              DisplacementFieldType;
	typedef typename DisplacementFieldType::ConstPointer    DisplacementFieldPointer;
	typedef typename DisplacementFieldType::PixelType       VectorType;
	//typedef typename VectorType::ValueType                  ValueType;
	typedef float ValueType;

	itkStaticConstMacro( Dimension, unsigned int, itkGetStaticConstMacro(VectorType::Dimension) );

	itkSetConstObjectMacro(Input,DisplacementFieldType);
	itkSetStringMacro(FileName);

	void Update() const {
		typedef itk::Image<ValueType,Dimension> FieldType;
		typename FieldType::Pointer out[Dimension];

		typename FieldType::SizeType outSize = m_Input->GetLargestPossibleRegion().GetSize();
		typename FieldType::SpacingType outSpacing = m_Input->GetSpacing();
		typename FieldType::DirectionType outDirection = m_Input->GetDirection();
		typename FieldType::PointType outOrigin = m_Input->GetOrigin();
		ValueType* buffer[3];

		for( size_t comp = 0; comp<Dimension; comp++ ){
			out[comp] = FieldType::New();
			out[comp]->SetRegions( outSize );
			out[comp]->SetSpacing( outSpacing );
			out[comp]->SetDirection( outDirection );
			out[comp]->SetOrigin( outOrigin );
			out[comp]->Allocate();
			out[comp]->FillBuffer(0.0);
			buffer[comp] = out[comp]->GetBufferPointer();
		}

		const VectorType* vectBuffer = m_Input->GetBufferPointer();
		size_t nVect = m_Input->GetLargestPossibleRegion().GetNumberOfPixels();
		for(size_t pix = 0; pix< nVect; pix++) {
			VectorType val = *(vectBuffer+pix);
			for( size_t comp=0; comp<Dimension; comp++) {
				*(buffer[comp]+pix) = static_cast<ValueType> (val[comp]);
			}
		}

		for( size_t comp=0; comp<Dimension; comp++) {
			std::stringstream ss;
			ss << m_FileName << "_cmp" << comp << ".nii.gz";
			typename itk::ImageFileWriter<FieldType>::Pointer w = itk::ImageFileWriter<FieldType>::New();
			w->SetInput( out[comp] );
			w->SetFileName( ss.str().c_str() );
			w->Update();
		}
	}

protected:
	DisplacementFieldComponentsFileWriter(){}
	~DisplacementFieldComponentsFileWriter(){}

	void PrintSelf( std::ostream& os, itk::Indent indent) const {
		Superclass::Print(os, indent);
	}

private:
	DisplacementFieldComponentsFileWriter( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	std::string m_FileName;
	DisplacementFieldPointer m_Input;
}; // End of class DisplacementFieldComponentsFileWriter
} // End of namespace


#endif /* DisplacementFieldComponentsFileWriter_H_ */
