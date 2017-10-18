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

#ifndef RSTKVTKPOLYDATAWRITER_H_
#define RSTKVTKPOLYDATAWRITER_H_

#include <itkVTKPolyDataWriter.h>

namespace rstk {
/** \class rstkVTKPolyDataWriter
 * \brief
 * Modifies the standard VTKPolyDataWriter behavior to write files
 * compatible with freeview (from FreeSurfer).
 * \ingroup ITKMesh
 */

template< class TInputMesh >
class VTKPolyDataWriter:
		public itk::VTKPolyDataWriter< TInputMesh >
{
public:
  /** Standard "Self" typedef. */
  typedef VTKPolyDataWriter                        Self;
  typedef itk::VTKPolyDataWriter< TInputMesh >     Superclass;
  typedef itk::SmartPointer< Self >                Pointer;
  typedef itk::SmartPointer< const Self >          ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VTKPolyDataWriter, itk::VTKPolyDataWriter);

  /** Hold on to the type information specified by the template parameters.
   */
  typedef typename Superclass::InputMeshType      InputMeshType;
  typedef typename Superclass::PixelType          PixelType;
  typedef typename Superclass::PointType          PointType;
  typedef typename Superclass::CellType           CellType;
  typedef typename Superclass::PointIdentifier    PointIdentifier;

  /** Some convenient typedefs. */
  typedef typename Superclass::ConstPointer InputMeshPointer;
  typedef typename Superclass::CellTraits   CellTraits;

  /** Define the triangular cell types which form the surface  */
  typedef typename Superclass::CellInterfaceType  CellInterfaceType;
  typedef typename Superclass::TriangleCellType   TriangleCellType;

  typedef typename Superclass::PointsContainer PointsContainer;
  typedef typename Superclass::CellsContainer  CellsContainer;

  typedef typename PointsContainer::ConstIterator PointIterator;
  typedef typename CellsContainer::ConstIterator  CellIterator;

  typedef typename CellType::PointIdIterator PointIdIterator;

protected:
  VTKPolyDataWriter(): Superclass(){}
  virtual ~VTKPolyDataWriter(){}

  virtual void GenerateData();
  void PrintSelf(std::ostream & os, itk::Indent indent) const {
	  Superclass::PrintSelf( os, indent );
  }

private:
  VTKPolyDataWriter(const Self &); //purposely not implemented
  void operator=(const Self &);    //purposely not implemented
};
}  // namespace rstk


#ifndef ITK_MANUAL_INSTANTIATION
#include "rstkVTKPolyDataWriter.hxx"
#endif


#endif /* RSTKVTKPOLYDATAWRITER_H_ */
