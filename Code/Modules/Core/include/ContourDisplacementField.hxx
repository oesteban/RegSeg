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

#ifndef CONTOURDISPLACEMENTFIELD_HXX_
#define CONTOURDISPLACEMENTFIELD_HXX_

#include "ContourDisplacementField.h"

namespace rstk {

template< typename TCoordPrecision, unsigned int VDimension>
ContourDisplacementField< TCoordPrecision, VDimension>::ContourDisplacementField() {
	//this->m_ReferenceMesh = InternalMeshType::New();

}

template< typename TCoordPrecision, unsigned int VDimension>
void
ContourDisplacementField< TCoordPrecision, VDimension>::
PrintSelf( std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf( os, indent );

	os << "* Internal Mesh: " << std::endl;
	//this->m_ReferenceMesh->Print( os, indent.GetNextIndent() );
}

template< typename TCoordPrecision, unsigned int VDimension>
void
ContourDisplacementField< TCoordPrecision, VDimension>::
Initialize() {
	typename Superclass::PointsContainerIterator d_it = this->GetPoints()->Begin();
	VectorType zeroVect = itk::NumericTraits<VectorType>::Zero;
	while (d_it != this->GetPoints()->End()) {
		this->SetPointData( d_it.Index(), zeroVect );
		++d_it;
	}
}

/*
template< typename TCoordPrecision, unsigned int VDimension>
void
ContourDisplacementField< TCoordPrecision, VDimension>::
SetReferenceMesh(InternalMeshType* mesh) {
	typename InternalMeshType::PointsContainerPointer points = mesh->GetPoints();
	typename InternalMeshType::PointsContainerIterator p_it =  points->Begin();

	typename InternalMeshType::PointDataContainerPointer container = mesh->GetPointData();
	typename InternalMeshType::PointDataContainerIterator d_it =     mesh->Begin();

	VectorType zeroVect();
	zeroVect.Fill( 0.0 );
	while (p_it != points->End()) {
		typename Superclass::PointIdentifier id = this->AddPoint( *p_it );
		this->SetPointData( id, zeroVect );
		++p_it;
	}

	typename InternalMeshType::CellsContainerPointer cells = mesh->GetCells();
	typename InternalMeshType::CellsContainerIterator c_it =  cells->Begin();

	while( c_it!= cells->End() ) {
		this->AddFace( )
		++c_it;
	}


}*/

}


#endif /* CONTOURDISPLACEMENTFIELD_HXX_ */
