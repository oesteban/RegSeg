/* --------------------------------------------------------------------------------------
 * File:    IDWMultivariateInterpolator.h
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
#include "IDWMultivariateInterpolator.txx"
#endif

#endif /* IDWMULTIVARIATEINTERPOLATOR_H_ */
