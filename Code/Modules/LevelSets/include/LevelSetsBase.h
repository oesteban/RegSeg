// --------------------------------------------------------------------------
// File:             LevelSetsBase.h
// Date:             27/10/2012
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
// This file is part of ACWERegistration-Debug@Debug
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

#ifndef LEVELSETSBASE_H_
#define LEVELSETSBASE_H_

#include <itkObject.h>
#include <itkNumericTraits.h>
#include <itkVector.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkQuadEdgeMesh.h>

namespace rstk {
/** \class LevelSetsBase
 *  \brief This class is a base for LevelSets definitions in the RSTK framework
 *
 *  This class is meant to provide an interface for LevelSets definitions. Given
 *  that RSTK is a registration framework, it provides the necessary elements to
 *  use the LevelSets function as a metric for non-linear registration.
 *
 *  LevelSetsBase requires the instantiation with a target image (DWI), a moving contour
 *  (surface) from T1 space and deformation fields definitions for fixed and moving
 *  objects. This implies to have a dense deformation field for DWI space and a sparse
 *  deformation field for the contour.
 *
 *  Derived classes must provide implementation for
 *  GetValue
 *  GetLevelSetsMap
 *
 *  \ingroup LevelSets
 *  \ingroup RSTK
 */


template <typename TReferenceImageType, typename TCoordRepType = double>
class LevelSetsBase: public itk::Object {
public:
	typedef LevelSetsBase                    Self;
	typedef itk::Object                      Superclass;
	typedef itk::SmartPointer<Self>          Pointer;
	typedef itk::SmartPointer< const Self >  ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro(LevelSetsBase, itk::Object);

	itkStaticConstMacro( Dimension, unsigned int, TReferenceImageType::ImageDimension );

	typedef double                                           MeasureType;
	typedef TCoordRepType                                    PointValueType;
	typedef itk::Vector< PointValueType, Dimension >         VectorType;

	typedef TReferenceImageType                              ReferenceImageType;
	typedef typename ReferenceImageType::PixelType           PixelType;
	typedef typename PixelType::ValueType                    PixelValueType;

	typedef itk::VectorLinearInterpolateImageFunction
			< ReferenceImageType >                           InterpolatorType;
	typedef typename InterpolatorType::Pointer               InterpolatorPointer;

	typedef itk::QuadEdgeMesh< VectorType, Dimension >       ContourDeformationType;
	typedef typename ContourDeformationType::PointType       PointType;
	typedef typename ContourDeformationType::Pointer         ContourDeformationPointer;
	typedef typename ContourDeformationType
			                     ::PointDataContainerPointer PointDataContainerPointer;

	typedef itk::Image< VectorType, Dimension >              DeformationFieldType;
	typedef typename DeformationFieldType::Pointer           DeformationFieldPointer;
	typedef typename DeformationFieldType::ConstPointer      DeformationFieldConstPointer;

	virtual MeasureType GetValue() const = 0;
	virtual void GetLevelSetsMap( DeformationFieldType* levelSetMap) const = 0;
protected:
	LevelSetsBase() { this->m_Value = itk::NumericTraits<MeasureType>::infinity(); }
	virtual ~LevelSetsBase() {}

	void PrintSelf(std::ostream & os, itk::Indent indent) const {
		Superclass::PrintSelf(os, indent);
		os << indent << "Value: " << m_Value << std::endl;
		os << std::endl;
	}

	mutable MeasureType m_Value;
	mutable DeformationFieldConstPointer m_DeformationField;
private:
	LevelSetsBase(const Self &);  //purposely not implemented
	void operator=(const Self &); //purposely not implemented

}; // end LevelSetsBase Class
} // end namespace rstk
#endif /* LEVELSETSBASE_H_ */
