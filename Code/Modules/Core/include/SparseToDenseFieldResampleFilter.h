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
template< class TInputMesh, class TDenseFieldType >
class SparseToDenseFieldResampleFilter: public MeshToImageFilter< TInputMesh, TDenseFieldType >
{
public:
	typedef SparseToDenseFieldResampleFilter                    Self;
	typedef TInputMesh                                          ContourType;
	typedef TDenseFieldType                                     FieldType;
	typedef MeshToImageFilter< ContourType, FieldType >         Superclass;
	typedef itk::SmartPointer< Self >                           Pointer;
	typedef itk::SmartPointer< const Self >                     ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( SparseToDenseFieldResampleFilter, itk::MeshToImageFilter );

	itkStaticConstMacro( Dimension, unsigned int, FieldType::ImageDimension );

	typedef typename FieldType::Pointer                    FieldPointer;
	typedef typename FieldType::ConstPointer               FieldConstPointer;
	typedef typename FieldType::PixelType                  VectorType;
	typedef typename FieldType::PointType                  FieldPointType;
	typedef typename FieldType::SpacingType                FieldSpacingType;
	typedef typename FieldType::IndexType                  FieldIndexType;
	typedef typename FieldType::DirectionType              FieldDirectionType;
	typedef typename FieldType::SizeType                   FieldSizeType;
	typedef typename VectorType::ValueType                 CoordinateValueType;
	typedef itk::ContinuousIndex<
		typename FieldPointType::ValueType, Dimension >    FieldContinuousIndexType;

	typedef typename itk::Point
			       < CoordinateValueType, Dimension>       PointType;

	typedef typename std::vector<PointType>                ControlPointList;

	typedef vnl_sparse_matrix< CoordinateValueType >       WeightsMatrix;
	typedef vnl_vector< CoordinateValueType >              GradientsVector;
	typedef vnl_matrix< PointType >                        DenseFieldMatrix;

	void AddControlPoints( const ContourType* prior );
	inline void AddControlPoint( const PointType p );

	itkGetMacro( ControlPoints, ControlPointList );


	/** Set the output image spacing. */
	itkSetMacro(FieldSpacing, FieldSpacingType);
	virtual void SetFieldSpacing(const double *values);

	/** Get the output image spacing. */
	itkGetConstReferenceMacro(FieldSpacing, FieldSpacingType);

	/** Set the output image origin. */
	itkSetMacro(FieldOrigin, FieldPointType);
	virtual void SetFieldOrigin(const double *values);

	/** Get the output image origin. */
	itkGetConstReferenceMacro(FieldOrigin, FieldPointType);

	itkSetMacro(FieldSize, FieldSizeType);
	itkGetConstReferenceMacro(FieldSize, FieldSizeType);

	/** Set the output direciton cosine matrix. */
	itkSetMacro(FieldDirection, FieldDirectionType);
	itkGetConstReferenceMacro(FieldDirection, FieldDirectionType);

	/** Set the start index of the output largest possible region.
	 * The default is an index of all zeros. */
	itkSetMacro(FieldStartIndex, FieldIndexType);

	/** Get the start index of the output largest possible region. */
	itkGetConstReferenceMacro(FieldStartIndex, FieldIndexType);

	virtual void CopyImageInformation( const FieldType* image );

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

	//InterpolatorPointerType                      m_Interpolator;
	FieldDirectionType                          m_FieldDirection;
	FieldPointType                              m_FieldOrigin;
	FieldSpacingType                            m_FieldSpacing;
	FieldIndexType                              m_FieldStartIndex;
	FieldSizeType                               m_FieldSize;
	FieldPointer								m_Field;

	ControlPointList                             m_ControlPoints;

	WeightsMatrix                                m_Phi;
	GradientsVector                              m_LevelSetVector[Dimension];
	size_t                                       m_N;
	size_t                                       m_k;
	bool                                         m_IsPhiInitialized;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "SparseToDenseFieldResampleFilter.hxx"
#endif


#endif /* SPARSETODENSEFIELDRESAMPLEFILTER_H_ */
