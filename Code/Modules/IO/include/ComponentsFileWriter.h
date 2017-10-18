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
