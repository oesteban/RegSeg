// --------------------------------------------------------------------------------------
// File:          WarpQEMeshFilter.h
// Date:          Jan 14, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

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
