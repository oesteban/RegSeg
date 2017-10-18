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

#ifndef IDWMULTIVARIATEINTERPOLATOR_H_
#define IDWMULTIVARIATEINTERPOLATOR_H_


#include "SparseMultivariateInterpolator.h"

using namespace itk;

namespace rstk {

/** \class IDWMultivariateInterpolator
 * \brief
 *
 * IDWMultivariateInterpolator is the base class for all process objects that output
 * image data and require a Mesh data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ImageFilters
 * \ingroup ITKMesh
 */
template< class TInputMesh, class TInputImage >
class IDWMultivariateInterpolator:public SparseMultivariateInterpolator< TInputMesh, TInputImage >
{
public:
	typedef IDWMultivariateInterpolator                    Self;
	typedef SparseMultivariateInterpolator
			                 < TInputMesh, TInputImage >   Superclass;
	typedef SmartPointer< Self >                           Pointer;
	typedef SmartPointer< const Self >                     ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( IDWMultivariateInterpolator, SparseMultivariateInterpolator );


	/** InputImageType typedef support. */
	typedef typename Superclass::InputImageType            InputImageType;
	typedef typename Superclass::PixelType                 PixelType;
	typedef typename Superclass::ValueType                 ValueType;
	typedef typename Superclass::RealType                  RealType;
	typedef typename Superclass::PointType                 PointType;
	typedef typename Superclass::IndexType                 IndexType;
	typedef typename Superclass::ContinuousIndexType       ContinuousIndexType;
	typedef typename Superclass::OutputType                OutputType;
	typedef typename Superclass::CoordRepType              CoordRepType;

	/** InputMeshType typedef support */
	typedef typename Superclass::InputMeshType             InputMeshType;
	typedef typename Superclass::InputMeshPointer          InputMeshPointer;
	typedef typename Superclass::InputMeshConstPointer     InputMeshConstPointer;


	itkStaticConstMacro( InputImageDimension, unsigned int, TInputImage::ImageDimension );

	  /** Interpolate the image at a continuous index position
	   *
	   * Returns the interpolated image intensity at a
	   * specified index position. No bounds checking is done.
	   * The point is assume to lie within the image buffer.
	   *
	   * Subclasses must override this method.
	   *
	   * ImageFunction::IsInsideBuffer() can be used to check bounds before
	   * calling the method. */
	 OutputType Evaluate( const PointType & point) const;


protected:
	IDWMultivariateInterpolator();
	~IDWMultivariateInterpolator(){};

	virtual void PrintSelf( std::ostream &os, Indent indent) const;

private:
	IDWMultivariateInterpolator( const Self& ); // purposely not implemented
	void operator=( const Self& );                 // purposely not implemented

	double m_pFactor;

};
}


#ifndef ITK_MANUAL_INSTANTIATION
#include "IDWMultivariateInterpolator.hxx"
#endif


#endif /* IDWMULTIVARIATEINTERPOLATOR_H_ */
