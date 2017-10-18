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

#ifndef SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_HXX_
#define SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_HXX_

#include "MultilabelBinarizeMeshFilter.h"
#include <itkProcessObject.h>
#include <itkImageRegionIterator.h>

namespace rstk
{
/** Constructor */
template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::MultilabelBinarizeMeshFilter():
 m_NumberOfMeshes(0) {
	this->SetNumberOfRequiredInputs(1);
	m_Size.Fill(0);
	m_Index.Fill(0);
	m_Spacing.Fill(0.0);
	m_Origin.Fill(0.0);
	m_Direction.GetVnlMatrix().set_identity();

	this->m_PreThreader = itk::MultiThreader::New();
	this->m_NumberOfThreads = this->m_PreThreader->GetNumberOfThreads();
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
void
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::SetOutputReference(const itk::ImageBase<Dimension>* reference) {
	this->SetSize(reference->GetLargestPossibleRegion().GetSize());
	this->SetOrigin(reference->GetOrigin());
	this->SetDirection(reference->GetDirection());
	this->SetSpacing(reference->GetSpacing());
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
void
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::GenerateOutputInformation() {
	// Get the output pointer
	OutputImagePointer output = this->GetOutput();
	m_NumberOfRegions = m_NumberOfMeshes + 1;

	typename OutputImageType::RegionType region;
	region.SetSize(m_Size);
	region.SetIndex(m_Index);

	output->SetLargestPossibleRegion(region); //
	output->SetBufferedRegion(region);        // set the region
	output->SetRequestedRegion(region);       //
	output->SetSpacing(m_Spacing);            // set spacing
	output->SetOrigin(m_Origin);              //   and origin
	output->SetDirection(m_Direction);        // direction cosines
	output->SetNumberOfComponentsPerPixel(m_NumberOfRegions);
	output->Allocate();


	OutputPixelType zero;
	zero.SetSize(m_NumberOfRegions);
	zero.Fill(0);
	zero[m_NumberOfMeshes] = 1;
	output->FillBuffer(zero);

	m_OutputSegmentation = OutputComponentType::New();
	m_OutputSegmentation->SetLargestPossibleRegion(region); //
	m_OutputSegmentation->SetBufferedRegion(region);        // set the region
	m_OutputSegmentation->SetRequestedRegion(region);       //
	m_OutputSegmentation->SetSpacing(m_Spacing);            // set spacing
	m_OutputSegmentation->SetOrigin(m_Origin);              //   and origin
	m_OutputSegmentation->SetDirection(m_Direction);        // direction cosines
	m_OutputSegmentation->Allocate();
	m_OutputSegmentation->FillBuffer(m_NumberOfMeshes);
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
void
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::BinarizeThreaded(size_t num) {
	this->m_Components[num]->Update();
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
void
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::SplitRequestedFilters(itk::ThreadIdType id, itk::ThreadIdType num, std::vector<size_t>& res) {
	size_t total = this->m_Components.size();

	if (num > total)
		itkExceptionMacro(<< "total number of threads overs number of filters")

	res.push_back(id);

	for(int rem = total; rem > num; ) {
		size_t nextId = id + num;
		if( nextId < total )
			res.push_back(nextId);
		else
			return;

		rem-=num;
	}
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
ITK_THREAD_RETURN_TYPE
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::BinarizeThreaderCallback(void *arg) {
	ThreadStruct *str;
	itk::ThreadIdType total, threadId, threadCount;

	threadId = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
	threadCount = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;
	str = (ThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);

	std::vector<size_t> filters;
	str->Filter->SplitRequestedFilters(threadId, threadCount, filters);

	for (size_t id = 0; id < filters.size(); id++)
		str->Filter->BinarizeThreaded(filters[id]);

	return ITK_THREAD_RETURN_VALUE;
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
void
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::BeforeThreadedGenerateData() {
	// Get the input and output pointers
	OutputImagePointer output = this->GetOutput();
	m_Components.clear();

	// find the actual number of threads
	m_NumberOfThreads = this->GetNumberOfThreads();
	if ( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 ) {
		m_NumberOfThreads = vnl_math_min( this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
	}

	long nbOfThreads = (m_NumberOfMeshes > m_NumberOfThreads)?m_NumberOfThreads:m_NumberOfMeshes;
	struct ThreadStruct str;
	str.Filter = this;

	for (size_t idx = 0; idx < m_NumberOfMeshes; idx++ ) {
		BinarizeMeshFilterPointer meshFilter = BinarizeMeshFilterType::New();
		meshFilter->SetSpacing( m_Spacing );
		meshFilter->SetDirection( m_Direction );
		meshFilter->SetOrigin( m_Origin );
		meshFilter->SetSize( m_Size );
		meshFilter->SetInput( GetInput(idx) );
		m_Components.push_back(meshFilter);
	}

	this->GetPreMultiThreader()->SetNumberOfThreads( nbOfThreads );
	this->GetPreMultiThreader()->SetSingleMethod( this->BinarizeThreaderCallback, &str );
	this->GetPreMultiThreader()->SingleMethodExecute();

	// number of threads can be constrained by the region size, so call the
	// SplitRequestedRegion
	// to get the real number of threads which will be used
	RegionType splitRegion;  // dummy region - just to call the following method
	nbOfThreads = this->SplitRequestedRegion(0, m_NumberOfThreads, splitRegion);
}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
void
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::ThreadedGenerateData(const RegionType & inputRegionForThread, itk::ThreadIdType threadId) {
	size_t nPix = inputRegionForThread.GetNumberOfPixels();
	itk::ProgressReporter progress( this, threadId, nPix );

	const OutputPixelValueType* compBuffer[m_NumberOfMeshes];
	for( size_t comp = 0; comp < m_NumberOfMeshes; comp++ ) {
		compBuffer[comp] = m_Components[comp]->GetOutput()->GetBufferPointer();
	}
	OutputComponentPointer ref = m_OutputSegmentation;
	OutputPixelValueType* segBuffer = ref->GetBufferPointer();

	itk::ImageRegionIterator< OutputImageType > outIt( this->GetOutput(), inputRegionForThread );
 	outIt.GoToBegin();

 	size_t pix;
 	OutputPixelValueType* outBuffer = this->GetOutput()->GetBufferPointer();
 	OutputPixelValueType v;
 	while ( !outIt.IsAtEnd() ) {
 		pix = ref->ComputeOffset(outIt.GetIndex());

		for( size_t comp = 0; comp < m_NumberOfMeshes; comp++ ) {
			v = *( compBuffer[comp] + pix );

			if( v > 0 ) {
				*( outBuffer + pix * m_NumberOfRegions + comp ) = 1;
				*( outBuffer + pix * m_NumberOfRegions + m_NumberOfMeshes ) = 0;
				*( segBuffer + pix ) = comp;
				break;
			}
		}
		progress.CompletedPixel();
		++outIt;
 	}

}

template< typename TInputMesh, typename TOutputPixelType, unsigned int VDimension >
typename MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >::InputMeshType *
MultilabelBinarizeMeshFilter< TInputMesh, TOutputPixelType, VDimension >
::GetInput(size_t idx) {
  return itkDynamicCastInDebugMode< TInputMesh * >
         ( this->ProcessObject::GetInput(idx) );
}

} // namespace rstk
#endif /* SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_MULTILABELBINARIZEMESHFILTER_HXX_ */
