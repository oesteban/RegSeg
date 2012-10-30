// --------------------------------------------------------------------------
// File:             SparseMultivariateInterpolator.h
// Date:             30/10/2012
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

#ifndef SPARSEMULTIVARIATEINTERPOLATOR_H_
#define SPARSEMULTIVARIATEINTERPOLATOR_H_


#include "itkImageFunction.h"
#include "itkVector.h"

/** \class SparseMultivariateInterpolator
 * \brief Base class for all sparse to dense vector field interpolators.
 *
 * SparseMultivariateInterpolator is the base for all MeshToImageFunctions that
 * interpolates image with vector pixel types. This function outputs
 * a return value of type Vector<double,Dimension>.
 *
 * This class is templated input image type and the coordinate
 * representation type.
 *
 * \warning This hierarchy of functions work only for images
 * with Vector-based pixel types. For scalar images use
 * InterpolateImageFunction.
 *
 * \sa InterpolateImageFunction
 * \ingroup ImageFunctions ImageInterpolators
 * \ingroup ITKImageFunction
 */
namespace rstk {

template<class TInputMesh, class TInputImage, class TCoordRep = double>
class SparseMultivariateInterpolator:
	public itk::ImageFunction<TInputImage,
	                          typename itk::NumericTraits<typename TInputImage::PixelType>::RealType,
	                          TCoordRep>
{
public:

	/** Standard class typedefs. */
	typedef SparseMultivariateInterpolator                            Self;
	typedef typename itk::NumericTraits
			<typename TInputImage::PixelType>::RealType               PixelRealType;
	typedef itk::ImageFunction<TInputImage, PixelRealType,TCoordRep>  Superclass;
	typedef itk::SmartPointer<Self>                                   Pointer;
	typedef itk::SmartPointer<const Self>                             ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro(SparseMultivariateInterpolator, ImageFunction);

	/** Extract the vector dimension from the pixel template parameter. */
	itkStaticConstMacro(VectorDimension, unsigned int, TInputImage::PixelType::Dimension);

	/** Dimension underlying input image. */
	itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);


	/** InputImageType typedef support. */
	typedef typename Superclass::InputImageType                  InputImageType;
	typedef typename InputImageType::PixelType                   PixelType;
	typedef typename PixelType::ValueType                        ValueType;
	typedef typename itk::NumericTraits<ValueType>::RealType     RealType;
	typedef typename Superclass::PointType                       PointType;
	typedef typename Superclass::IndexType                       IndexType;
	typedef typename Superclass::ContinuousIndexType             ContinuousIndexType;
	typedef typename Superclass::OutputType                      OutputType;
	typedef TCoordRep                                            CoordRepType;

	/** InputMeshType typedef support */
	typedef          TInputMesh                                  InputMeshType;
	typedef typename InputMeshType::Pointer                      InputMeshPointer;
	typedef typename InputMeshType::ConstPointer                 InputMeshConstPointer;
	typedef typename InputMeshType::PointIdentifier              InputMeshPointIdentifier;



	  /** Set the input mesh.
	   */
	  virtual void SetInputMesh(const InputMeshType *ptr);
	  /** Get the input image. */
	  const InputMeshType * GetInputMesh() const
	  { return m_Mesh.GetPointer(); }

	/** Returns the interpolated image intensity at a
	 * specified point position. No bounds checking is done.
	 * The point is assume to lie within the image buffer.
	 * ImageFunction::IsInsideBuffer() can be used to check bounds before
	 * calling the method. */
	virtual OutputType Evaluate(const PointType & point) const = 0;

	/** Interpolate the image at a continuous index position
	 *
	 * Returns the interpolated image intensity at a
	 * specified index position. No bounds checking is done.
	 * The point is assume to lie within the image buffer.
	 *
	 * Subclasses must override this method.
	 */
	virtual OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index) const {
		PointType point;

		this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(index, point);
		return (this->Evaluate(point));
	}

	/** Interpolate the image at an index position.
	 */
	virtual OutputType EvaluateAtIndex(const IndexType & index) const {
		PointType point;

		this->GetInputImage()->TransformIndexToPhysicalPoint(index, point);
		return (this->Evaluate(point));
	}

protected:
	SparseMultivariateInterpolator();
	~SparseMultivariateInterpolator(){}


	void PrintSelf(std::ostream & os, itk::Indent indent) const;

	InputMeshConstPointer     m_Mesh;
	InputMeshPointIdentifier  m_NumberOfMeshPoints;
	size_t                    m_NumberOfImagePixels;
private:
	SparseMultivariateInterpolator(const Self &); //purposely not implemented
	void operator=(const Self &); //purposely not implemented
};

}

#if ITK_TEMPLATE_TXX
#include "SparseMultivariateInterpolator.hxx"
#endif


#endif /* SPARSEMULTIVARIATEINTERPOLATOR_H_ */
