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
