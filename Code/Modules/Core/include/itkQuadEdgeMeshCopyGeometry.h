// --------------------------------------------------------------------------------------
// File:          itkQuadEdgeMeshCopyGeometry.h
// Date:          Jan 9, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of RegSeg
//
// RegSeg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RegSeg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RegSeg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef ITKQUADEDGEMESHCOPYGEOMETRY_H_
#define ITKQUADEDGEMESHCOPYGEOMETRY_H_

#include "itkMeshToMeshFilter.h"

namespace itk
{
/** \class itkQuadEdgeMeshCopyGeometry
 *  \brief Duplicates the content of a Mesh
 *
 * \author Alexandre Gouaillard, Leonardo Florez-Valencia, Eric Boix
 *
 * This implementation was contributed as a paper to the Insight Journal
 * http://hdl.handle.net/1926/306
 *
 * \ingroup ITKQuadEdgeMesh
 */
template< typename TInputMesh, typename TOutputMesh >
class ITK_EXPORT QuadEdgeMeshCopyGeometry:
  public MeshToMeshFilter< TInputMesh, TOutputMesh >
{
public:
  /** Basic types. */
  typedef QuadEdgeMeshCopyGeometry            Self;
  typedef MeshToMeshFilter< TInputMesh, TOutputMesh > Superclass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** Input types. */
  typedef TInputMesh                              InputMeshType;
  typedef typename InputMeshType::Pointer         InputMeshPointer;
  typedef typename InputMeshType::ConstPointer    InputMeshConstPointer;
  typedef typename InputMeshType::CoordRepType    InputCoordRepType;
  typedef typename InputMeshType::PointType       InputPointType;
  typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
  typedef typename InputMeshType::QEPrimal        InputQEPrimal;
  typedef typename InputMeshType::VectorType      InputVectorType;

  typedef typename InputMeshType::PointDataContainer InputPointDataContainer;
  typedef typename InputMeshType::CellDataContainer  InputCellDataContainer;

  typedef typename InputPointDataContainer::ConstPointer
  InputPointDataContainerConstPointer;
  typedef typename InputMeshType::PointsContainerConstIterator
  InputPointsContainerConstIterator;
  typedef typename InputMeshType::PointsContainerConstPointer
  InputPointsContainerConstPointer;
  typedef typename InputMeshType::CellsContainerConstIterator
  InputCellsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstPointer
  InputCellsContainerConstPointer;

  typedef typename InputMeshType::EdgeCellType    InputEdgeCellType;
  typedef typename InputMeshType::PolygonCellType InputPolygonCellType;
  typedef typename InputMeshType::PointIdList     InputPointIdList;
  typedef typename InputMeshType::CellTraits      InputCellTraits;
  typedef typename InputCellTraits::PointIdInternalIterator
  InputPointsIdInternalIterator;

  typedef typename InputQEPrimal::IteratorGeom InputQEIterator;

  /** Output types. */
  typedef TOutputMesh                              OutputMeshType;
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer    OutputMeshConstPointer;
  typedef typename OutputMeshType::CoordRepType    OutputCoordRepType;
  typedef typename OutputMeshType::PointType       OutputPointType;
  typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
  typedef typename OutputMeshType::QEPrimal        OutputQEPrimal;
  typedef typename OutputMeshType::VectorType      OutputVectorType;
  typedef typename OutputQEPrimal::IteratorGeom    OutputQEIterator;
  typedef typename OutputMeshType::PointsContainerIterator
  OutputPointsContainerIterator;
  typedef typename OutputMeshType::PointsContainerPointer
  OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerConstPointer
  OutputPointsContainerConstPointer;

  typedef typename OutputMeshType::PointDataContainer OutputPointDataContainer;
  typedef typename OutputMeshType::CellDataContainer  OutputCellDataContainer;

public:
  itkNewMacro(Self);
  itkTypeMacro(QuadEdgeMeshCopyGeometry, MeshToMeshFilter);

protected:
  QuadEdgeMeshCopyGeometry();
  virtual ~QuadEdgeMeshCopyGeometry() {}
  virtual void CopyMesh();
  virtual void CopyMeshGeometry();
  virtual void CopyMeshPoints();
  virtual void CopyMeshCells();
  virtual void CopyMeshEdgeCells();

private:
  QuadEdgeMeshCopyGeometry(const Self &); // Not impl.
  void operator=(const Self &);                   // Not impl.

  void Copy          (const TInputMesh *in, TOutputMesh *out);
  void CopyCells     (const TInputMesh *in, TOutputMesh *out);
  void CopyEdgeCells (const TInputMesh *in, TOutputMesh *out);
  void CopyPoints    (const TInputMesh *in, TOutputMesh *out);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshCopyGeometry.hxx"
#endif


#endif /* ITKQUADEDGEMESHCOPYGEOMETRY_H_ */
