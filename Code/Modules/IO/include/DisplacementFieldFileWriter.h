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
