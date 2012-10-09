/* --------------------------------------------------------------------------------------
 * File:    SparseDeformationField.h
 * Date:    Mar 27, 2012
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

#ifndef SPARSEDEFORMATIONFIELD_H_
#define SPARSEDEFORMATIONFIELD_H_

#include <vector>
#include <itkObject.h>
#include <itkVector.h>
#include <itkPointSet.h>

#include <vector>
#include <itkObject.h>
#include <itkVector.h>
#include <itkPointSet.h>

namespace rstk {

template <class TReferenceSurface, class TCoordinatesPrecision = float >
class SparseDeformationField: public itk::Object
{
public:
	/** Standard class typedefs */
	typedef SparseDeformationField                     Self;
	typedef itk::Object                                Superclass;
	typedef itk::SmartPointer<Self>                    Pointer;
	typedef itk::SmartPointer<const Self>              ConstPointer;

	typedef TCoordinatesPrecision                      PrecisionType;

	itkStaticConstMacro( Dimension, unsigned int,  TReferenceSurface::PointDimension );

	itkTypeMacro(SparseDeformationField, Object);
	itkNewMacro( Self );

	typedef itk::Vector< PrecisionType, Dimension >                ParameterType;
	typedef itk::PointSet< ParameterType, Dimension >              SparseVectorFieldType;
	typedef typename SparseVectorFieldType::Pointer                SparseVectorFieldPointer;
	typedef typename SparseVectorFieldType::ConstPointer           SparseVectorFieldConstPointer;
	typedef typename SparseVectorFieldType::PointsContainer        ParameterContainer;
	typedef typename SparseVectorFieldType::PointIdentifier        ParameterIdentifier; // signed long

	typedef TReferenceSurface                                      ReferenceSurfaceType;
	typedef typename ReferenceSurfaceType::Pointer                 ReferenceSurfacePointer;
	typedef typename ReferenceSurfaceType::ConstPointer            ReferenceSurfaceConstPointer;
	typedef typename ReferenceSurfaceType::CellType                CellType;
	typedef typename ReferenceSurfaceType::CellsContainer          CellsContainer;
	typedef typename ReferenceSurfaceType::CellsContainerConstIterator CellsContainerConstIterator;





	inline ParameterType GetK( ParameterIdentifier k ) const {
		ParameterType p;
		m_Field->GetPointData( k, &p );
		return p;
		//return *( m_Field->GetPointData() + k);
	}

	inline void SetK( const ParameterIdentifier k, const ParameterType& p) {
		m_Field->SetPointData(k, p);
		//*( m_Field->GetPointData() + k) = p;
		//for ( ParameterIdentifier i = 0; i < m_ParameterSize; i++)
		//	*(this->m_Field->GetPointData() + k + i) = p[i];
	}

	void SetReferenceSurface( const ReferenceSurfaceType * ref );
	itkGetConstObjectMacro( ReferenceSurface, ReferenceSurfaceType );

	const ReferenceSurfaceType*  GetTransformedSurface( );

protected:
	SparseDeformationField();
	~SparseDeformationField(){};

private:
	SparseDeformationField( const Self& );    //purposely not implemented
	void operator=(const Self& );            //purposely not implemented

	SparseVectorFieldPointer                 m_Field;
	ReferenceSurfaceConstPointer             m_ReferenceSurface;
	ReferenceSurfacePointer                  m_TransformedPointSet;
	ParameterIdentifier                      m_NumberOfParameters;
};

}

#include "SparseDeformationField.txx"

#endif /* SPARSEDEFORMATIONFIELD_H_ */
