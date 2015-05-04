// --------------------------------------------------------------------------------------
// File:          MultilabelBinarizeMeshFilter.hxx
// Date:          Jan 9, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2015, code@oscaresteban.es (Oscar Esteban)
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
