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
