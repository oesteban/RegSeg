// --------------------------------------------------------------------------------------
// File:          MultilabelBinarizeMeshFilter.h
// Date:          Jan 9, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2015, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWEReg
//
// ACWEReg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWEReg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWEReg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_H_
#define SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_H_

#include <vector>
#include <itkImageBase.h>
#include <itkImageSource.h>
#include <itkVectorImage.h>
#include <itkTriangleMeshToBinaryImageFilter.h>

namespace rstk {
template< typename TInputMesh, typename TOutputPixelType = unsigned char, unsigned int VDimension = 3 >
class MultilabelBinarizeMeshFilter: public itk::ImageSource< itk::VectorImage< TOutputPixelType, VDimension > >
{
public:
	  /** Standard class typedefs. */
	  typedef MultilabelBinarizeMeshFilter Self;
	  typedef TOutputPixelType                                    OutputPixelValueType;
	  typedef itk::VectorImage< TOutputPixelType, VDimension >    OutputImageType;

	  typedef itk::ImageSource< OutputImageType >                 Superclass;
	  typedef itk::SmartPointer< Self >                           Pointer;
	  typedef itk::SmartPointer< const Self >                     ConstPointer;

	  /** Method for creation through the object factory. */
	  itkNewMacro(Self);

	  /** Run-time type information (and related methods). */
	  itkTypeMacro(MultilabelBinarizeMeshFilter, itk::ImageSource);

	  itkStaticConstMacro( Dimension, unsigned int, VDimension );

	  typedef typename OutputImageType::Pointer                   OutputImagePointer;
	  typedef typename OutputImageType::PixelType                 OutputPixelType;
	  typedef typename OutputImageType::PointType                 PointType;
	  typedef typename OutputImageType::IndexType                 IndexType;
	  typedef typename OutputImageType::SizeType                  SizeType;
	  typedef typename OutputImageType::RegionType                RegionType;
	  typedef typename OutputImageType::ValueType                 ValueType;
	  typedef typename OutputImageType::SpacingType               SpacingType;
	  typedef typename OutputImageType::DirectionType             DirectionType;
	  typedef typename Superclass::OutputImageRegionType          OutputImageRegionType;

	  /** Some convenient typedefs. */
	  typedef TInputMesh                                          InputMeshType;
	  typedef typename InputMeshType::Pointer                     InputMeshPointer;
	  typedef typename InputMeshType::ConstPointer                InputMeshConstPointer;
	  typedef typename InputMeshType::PointType                   InputPointType;
	  typedef typename InputMeshType::PixelType                   InputPixelType;
	  typedef typename InputMeshType::MeshTraits::CellTraits      InputCellTraitsType;
	  typedef typename InputMeshType::CellType                    CellType;
	  typedef typename InputMeshType::CellsContainerPointer       CellsContainerPointer;
	  typedef typename InputMeshType::CellsContainerIterator      CellsContainerIterator;
	  typedef typename std::vector<InputMeshPointer>              InputMeshContainer;

	  typedef typename InputMeshType::PointsContainer             InputPointsContainer;
	  typedef typename InputPointsContainer::Pointer              InputPointsContainerPointer;
	  typedef typename InputPointsContainer::Iterator             InputPointsContainerIterator;

	  typedef itk::Image< OutputPixelValueType, VDimension >      OutputComponentType;
	  typedef typename OutputComponentType::Pointer               OutputComponentPointer;
	  typedef itk::TriangleMeshToBinaryImageFilter
				          <InputMeshType, OutputComponentType>	  BinarizeMeshFilterType;
	  typedef typename BinarizeMeshFilterType::Pointer            BinarizeMeshFilterPointer;


	  // Set/Get spacing
	  itkSetMacro(Spacing, SpacingType);
	  itkGetConstReferenceMacro(Spacing, SpacingType);

	  // Set/Get Direction
	  itkSetMacro(Direction, DirectionType);
	  itkGetConstMacro(Direction, DirectionType);

	  // Set/get origin
	  itkSetMacro(Origin, PointType);
	  itkGetConstReferenceMacro(Origin, PointType);

	  /** Set/Get Index */
	  itkSetMacro(Index, IndexType);
	  itkGetConstMacro(Index, IndexType);

	  /** Set/Get Size */
	  itkSetMacro(Size, SizeType);
	  itkGetConstMacro(Size, SizeType);

	  void SetOutputReference(const itk::ImageBase<Dimension>* reference);

	  InputMeshType * GetInput(size_t idx);
	  void SetInputs(const InputMeshContainer & cont) {
		  for(size_t c = 0; c < cont.size(); c++) {
			  InputMeshPointer ptr = cont[c];
			  this->ProcessObject::SetNthInput(c, ptr);
			  m_NumberOfMeshes++;
		  }
	  }

	  itkGetObjectMacro(OutputSegmentation, OutputComponentType)
	  itkGetConstObjectMacro(OutputSegmentation, OutputComponentType)
protected:
	  void BeforeThreadedGenerateData();
	  void ThreadedGenerateData(const RegionType & inputRegionForThread, itk::ThreadIdType threadId);
	  virtual void GenerateOutputInformation();

	  MultilabelBinarizeMeshFilter();
	  ~MultilabelBinarizeMeshFilter() {}
	  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const { Superclass::PrintSelf(os, indent); }

private:
	  MultilabelBinarizeMeshFilter(const Self &); //purposely not implemented
	  void operator=(const Self &);                  //purposely not implemented

	  IndexType m_Index;
	  SizeType m_Size;
	  SpacingType m_Spacing;
	  PointType m_Origin;
	  DirectionType m_Direction;

	  size_t m_NumberOfMeshes;
	  size_t m_NumberOfRegions;
	  std::vector< OutputComponentPointer > m_Components;
	  OutputComponentPointer m_OutputSegmentation;
}; // class

} // namespace rstk


#ifndef ITK_MANUAL_INSTANTIATION
#include "MultilabelBinarizeMeshFilter.hxx"
#endif

#endif /* SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_H_ */
