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
