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

	void PrintSelf( std::ostream& os, itk::Indent indent) const override {
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
