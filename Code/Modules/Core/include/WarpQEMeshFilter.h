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

#ifndef WARPQEMESHFILTER_H_
#define WARPQEMESHFILTER_H_

#include <itkMeshToMeshFilter.h>

namespace rstk {
/** \class WarpQEMeshFilter
 * \brief
 *
 * WarpQEMeshFilter applies a sparse deformation field (vector field defined at mesh point
 * locations).
 *
 * Locations of the sparse field are serialized vectors, one per dimension.
 *
 */

template< class TInputMesh, class TOutputMesh, class TVectorField >
class WarpQEMeshFilter: public itk::MeshToMeshFilter< TInputMesh, TOutputMesh >
{
public:
	/** Standard class typedefs. */
	typedef WarpQEMeshFilter                                       Self;
	typedef itk::MeshToMeshFilter< TInputMesh, TOutputMesh >       Superclass;
	typedef itk::SmartPointer< Self >                              Pointer;
	typedef itk::SmartPointer< const Self >                        ConstPointer;

	typedef TInputMesh                                             InputMeshType;
	typedef typename InputMeshType::Pointer                        InputMeshPointer;

	typedef TOutputMesh                                            OutputMeshType;
	typedef typename OutputMeshType::Pointer                       OutputMeshPointer;

	/** Type for representing coordinates */
	typedef typename TInputMesh::CoordRepType                      CoordRepType;

	/** Deformation field typedef support. */
	typedef TVectorField                                           FieldType;
	typedef typename FieldType::Pointer                            FieldPointer;
	typedef typename FieldType::ConstPointer                       FieldConstPointer;
	typedef typename FieldType::VectorType                         FieldVectorType;

	/** Method for creating object using the factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( WarpQEMeshFilter, itk::MeshToMeshFilter );

	/** Deformation field, get/set. */
	void SetField( const FieldType *field );
	const FieldType* GetField() const;

protected:
	WarpQEMeshFilter();
	~WarpQEMeshFilter() {}
	void PrintSelf( std::ostream & os, itk::Indent indent ) const;

	/** Generate Requested data */
	virtual void GenerateData();

private:
	WarpQEMeshFilter( const WarpQEMeshFilter & ); // purposely not implemented
	void operator=( const WarpQEMeshFilter & );   // purposely not implemented
};
} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "WarpQEMeshFilter.hxx"
#endif


#endif /* WARPQEMESHFILTER_H_ */
