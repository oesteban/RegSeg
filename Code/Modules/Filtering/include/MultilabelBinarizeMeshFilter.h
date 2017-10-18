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

	  typedef itk::ProcessObject                                  ProcessObject;


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

	  struct ThreadStruct { Pointer Filter; };
	  itk::MultiThreader * GetPreMultiThreader() const { return m_PreThreader; }
protected:
	  void BeforeThreadedGenerateData();
	  void ThreadedGenerateData(const RegionType & inputRegionForThread, itk::ThreadIdType threadId);
	  virtual void GenerateOutputInformation();

	  static ITK_THREAD_RETURN_TYPE BinarizeThreaderCallback(void *arg);
	  void BinarizeThreaded(size_t num);
	  void SplitRequestedFilters(itk::ThreadIdType id, itk::ThreadIdType num, std::vector<size_t>& res);

	  MultilabelBinarizeMeshFilter();
	  ~MultilabelBinarizeMeshFilter() {}
	  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const { Superclass::PrintSelf(os, indent); }

		/** Support processing data in multiple threads. */
		itk::MultiThreader::Pointer m_PreThreader;
		itk::ThreadIdType           m_NumberOfThreads;
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
	  std::vector< BinarizeMeshFilterPointer > m_Components;
	  OutputComponentPointer m_OutputSegmentation;
}; // class

} // namespace rstk


#ifndef ITK_MANUAL_INSTANTIATION
#include "MultilabelBinarizeMeshFilter.hxx"
#endif

#endif /* SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_H_ */
