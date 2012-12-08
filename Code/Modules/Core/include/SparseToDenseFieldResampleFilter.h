// --------------------------------------------------------------------------
// File:             SparseToDenseFieldResampleFilter.h
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

#ifndef SPARSETODENSEFIELDRESAMPLEFILTER_H_
#define SPARSETODENSEFIELDRESAMPLEFILTER_H_

#include "MeshToImageFilter.h"
#include "SparseMultivariateInterpolator.h"
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>


namespace rstk {

/** \class SparseToDenseFieldResampleFilter
 * \brief
 *
 * SparseToDenseFieldResampleFilter produces a dense deformation field interpolated
 * from a sparse (mesh) deformation field.
 *
 * \ingroup ImageFilters
 * \ingroup ITKMesh
 */
template< class TInputMesh, class TOutputImage >
class SparseToDenseFieldResampleFilter: public MeshToImageFilter< TInputMesh, TOutputImage >
{
public:
	typedef SparseToDenseFieldResampleFilter                    Self;
	typedef MeshToImageFilter< TInputMesh, TOutputImage >       Superclass;
	typedef itk::SmartPointer< Self >                           Pointer;
	typedef itk::SmartPointer< const Self >                     ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( SparseToDenseFieldResampleFilter, itk::MeshToImageFilter );

	typedef typename Superclass::InputMeshType             InputMeshType;
	typedef typename Superclass::InputMeshPointer          InputMeshPointer;
	typedef typename Superclass::InputMeshConstPointer     InputMeshConstPointer;
	typedef typename Superclass::OutputImageType           OutputImageType;
	typedef typename Superclass::OutputImagePointer        OutputImagePointer;
	typedef typename Superclass::OutputImageConstPointer   OutputImageConstPointer;
	typedef typename Superclass::OutputImageSizeType       OutputImageSizeType;
	typedef typename Superclass::ValueType                 ValueType;

	typedef typename OutputImageType::PixelType            OutputPixelType;
	typedef typename OutputImageType::PointType            OutputPointType;
	typedef typename OutputImageType::SpacingType          OutputSpacingType;
	typedef typename OutputImageType::IndexType            OutputIndexType;
	typedef typename OutputImageType::DirectionType        OutputDirectionType;
	typedef typename OutputImageType::RegionType           OutputRegionType;
	typedef typename OutputPixelType::ValueType            OutputVectorValueType;

	typedef vnl_sparse_matrix< double >                    WeightsMatrix;
	typedef vnl_vector< double >                           SpeedsVector;
	typedef vnl_matrix< OutputPixelType >                  DenseFieldMatrix;

	typedef SparseMultivariateInterpolator
			< TInputMesh, TOutputImage >                   InterpolatorType;
	typedef typename InterpolatorType::Pointer             InterpolatorPointerType;


	itkStaticConstMacro( OutputImageDimension, unsigned int, TOutputImage::ImageDimension );


	virtual void SetInput( const InputMeshType *input ) {
		this->SetInput( 0, input );
	}
	virtual void SetInput( unsigned int, const InputMeshType *input) {
		this->m_Mesh = input;
	}

	const InputMeshType*   GetInput(void) { return this->GetInput(0); }
	const InputMeshType*   GetInput(unsigned int idx) { return this->m_Mesh; }
	//const OutputImageType* GetOutput( void ) { return this->GetOutput(0); }

	itkSetConstObjectMacro( ShapePrior, InputMeshType );
	itkGetConstObjectMacro( ShapePrior, InputMeshType );

	/** Set the interpolator function.  The default is
	 * itk::VectorLinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>.  */
	itkSetObjectMacro(Interpolator, InterpolatorType);

	/** Get a pointer to the interpolator function. */
	itkGetConstObjectMacro(Interpolator, InterpolatorType);

	/** Set the size of the output image. */
	itkSetMacro(OutputSize, OutputImageSizeType);

	/** Get the size of the output image. */
	itkGetConstReferenceMacro(OutputSize, OutputImageSizeType);

	/** Set the pixel value when a transformed pixel is outside of the
	 * image.  The default default pixel value is 0. */
	itkSetMacro(DefaultPixelValue, OutputPixelType);

	/** Get the pixel value when a transformed pixel is outside of the image */
	itkGetConstMacro(DefaultPixelValue, OutputPixelType);

	/** Set the output image spacing. */
	itkSetMacro(OutputSpacing, OutputSpacingType);
	virtual void SetOutputSpacing(const double *values);

	/** Get the output image spacing. */
	itkGetConstReferenceMacro(OutputSpacing, OutputSpacingType);

	/** Set the output image origin. */
	itkSetMacro(OutputOrigin, OutputPointType);
	virtual void SetOutputOrigin(const double *values);

	/** Get the output image origin. */
	itkGetConstReferenceMacro(OutputOrigin, OutputPointType);

	/** Set the output direciton cosine matrix. */
	itkSetMacro(OutputDirection, OutputDirectionType);
	itkGetConstReferenceMacro(OutputDirection, OutputDirectionType);

	/** Set the start index of the output largest possible region.
	 * The default is an index of all zeros. */
	itkSetMacro(OutputStartIndex, OutputIndexType);

	/** Get the start index of the output largest possible region. */
	itkGetConstReferenceMacro(OutputStartIndex, OutputIndexType);

	virtual void CopyImageInformation( const OutputImageType* image );

	/** VectorResampleImageFilter produces an image which is a different size
	 * than its input.  As such, it needs to provide an implementation
	 * for GenerateOutputInformation() in order to inform the pipeline
	 * execution model.  The original documentation of this method is
	 * below. \sa ProcessObject::GenerateOutputInformaton() */
	virtual void GenerateOutputInformation();
	virtual void UpdateOutputInformation();

	/** Method Compute the Modified Time based on changed to the components. */
	size_t GetMTime(void) const;

protected:
	SparseToDenseFieldResampleFilter();
	~SparseToDenseFieldResampleFilter(){};

	virtual void GenerateData();

	virtual void PrintSelf( std::ostream &os, Indent indent) const;
private:
	SparseToDenseFieldResampleFilter( const Self& ); // purposely not implemented
	void operator=( const Self& );                 // purposely not implemented

	ValueType                       m_DefaultPixelValue;
	InterpolatorPointerType         m_Interpolator;
	OutputDirectionType             m_OutputDirection;
	OutputPointType                 m_OutputOrigin;
	OutputSpacingType               m_OutputSpacing;
	OutputIndexType                 m_OutputStartIndex;
	OutputImageSizeType             m_OutputSize;

	InputMeshConstPointer           m_ShapePrior;

	WeightsMatrix                   m_Phi;
	SpeedsVector                    m_LevelSetVector[OutputImageDimension];
	SpeedsVector                    m_OutputMatrix;
	size_t                          m_N;
	size_t                          m_k;
	bool                            m_IsPhiInitialized;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "SparseToDenseFieldResampleFilter.hxx"
#endif


#endif /* SPARSETODENSEFIELDRESAMPLEFILTER_H_ */
