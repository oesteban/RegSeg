// --------------------------------------------------------------------------------------
// File:          rstkVTKPolyDataWriter.h
// Date:          Mar 17, 2014
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
