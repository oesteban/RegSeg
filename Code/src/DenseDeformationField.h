/* --------------------------------------------------------------------------------------
 * File:    DenseDeformationField.h
 * Date:    Mar 27, 2012
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

#ifndef DENSEDEFORMATIONFIELD_H_
#define DENSEDEFORMATIONFIELD_H_

#include <vector>
#include <itkObject.h>
#include <itkVector.h>
#include <itkImage.h>

template <class TReferenceImage, class TRealType = float >
class DenseDeformationField: public itk::Object
{
public:
	/** Standard class typedefs */
	typedef DenseDeformationField                      Self;
	typedef itk::Object                                Superclass;
	typedef SmartPointer<Self>                         Pointer;
	typedef SmartPointer<const Self>                   ConstPointer;

	typedef TRealType                                  ScalarType;

	itkTypeMacro(DenseDeformationField, Object);
	itkStaticConstMacro( Dimension, unsigned int,  TReferenceImage::ImageDimension );

	typedef itk::Vector< ScalarType, Dimension >       ParameterType;
	typedef itk::Image< ParameterType, Dimension >     FieldType;
	typedef typename FieldType::OffsetValueType        ParameterIdentifier; // signed long
	typedef typename FieldType::OffsetType             ParameterOffsetType;

	inline ParameterType GetK( ParameterIdentifier k ) {
		return *(this->m_Field->GetBufferPointer() + k);
	}

	inline void SetK( ParameterIdentifier k, ParameterType& p) {
		for ( ParameterIdentifier i = 0; i < m_ParameterSize; i++)
			*(this->m_Field->GetBufferPointer() + k + i) = p[i];
	}


protected:
	DenseDeformationField();
	~DenseDeformationField(){};

private:
	DenseDeformationField( const Self& );    //purposely not implemented
	void operator=(const Self& );            //purposely not implemented

	FieldType                                m_Field;
	ParameterIdentifier                      m_NumberOfParameters;
	size_t                                   m_ParameterSize;

};


#endif /* DENSEDEFORMATIONFIELD_H_ */
