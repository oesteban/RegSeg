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

#ifndef CONTOURDISPLACEMENTFIELD_H_
#define CONTOURDISPLACEMENTFIELD_H_

// Include headers
#include <itkObject.h>
#include <itkVector.h>
#include <itkMesh.h>
#include <itkQuadEdgeMesh.h>

// Namespace declaration

namespace rstk {
/** \class ContourDisplacementField
 *  \brief This class implements sparse displacement field defined along a reference mesh (discrete 2D-manifold)
 *
 *  ContourDisplacementField implements a QuadEdgeMesh in order to having available all the QE filters (specifically,
 *  the Normals filter). It is built upon a FreeSurfer's mesh set via SetReferenceMesh(). ContourDisplacementField
 *  process 3D deformation surfaces and fields by default.
 *
 *  \ingroup RSTKCore
 */

template< typename TCoordPrecision = double, unsigned int VDimension = 3u >
class ContourDisplacementField: public itk::QuadEdgeMesh< itk::Vector< TCoordPrecision, VDimension >, VDimension > {
public:
	typedef ContourDisplacementField      Self;
	typedef itk::QuadEdgeMesh< itk::Vector< TCoordPrecision, VDimension >, VDimension >
	                                      Superclass;
	typedef itk::SmartPointer<Self>       Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro( ContourDisplacementField, itk::Object );
	itkNewMacro( Self );

	typedef typename Superclass::PixelType                      VectorType;
	typedef itk::Mesh<TCoordPrecision, VDimension >             InternalMeshType;
	typedef typename InternalMeshType::Pointer                  InternalMeshPointer;

	//itkSetObjectMacro( ReferenceMesh, InternalMeshType );
	//itkGetConstObjectMacro( ReferenceMesh, InternalMeshType );

	//void SetReferenceMesh( InternalMeshType* mesh );
	void Initialize();

protected:
	ContourDisplacementField();
	~ContourDisplacementField() {}

	void PrintSelf( std::ostream& os, itk::Indent indent) const;

private:
	ContourDisplacementField( const Self &); // purposely not implemented
	void operator=(const Self &); // purposely not implemented

	//InternalMeshPointer m_ReferenceMesh;

}; // End of class ContourDisplacementField
} // End of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "ContourDisplacementField.hxx"
#endif

#endif /* CONTOURDISPLACEMENTFIELD_H_ */
