// --------------------------------------------------------------------------
// File:             LevelSetsGenerateTestObjects.cxx
// Date:             31/10/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWEReg
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "./"
// data source http://code.google.com/p/v3dsolver/source/browse/data/
#endif

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkBinaryMask3DMeshSource.h>
#include <itkAddImageFilter.h>
#include "MahalanobisLevelSets.h"

using namespace rstk;


const unsigned int Dimension = 3;

int main(int argc, char *argv[]) {
	typedef float                                PixelType;
	typedef itk::Vector<float, 2u>               VectorPixelType;
	typedef itk::Image<PixelType, Dimension>            ImageType;
	typedef MahalanobisLevelSets<itk::Image<VectorPixelType, Dimension> >      LevelSetsType;

	typedef typename LevelSetsType::ContourDeformationType     ContourDeformationType;
	typedef typename ContourDeformationType::Pointer           ContourDeformationPointer;
	typedef typename LevelSetsType::VectorType                 VectorType;

	typedef itk::BinaryMask3DMeshSource< itk::Image<PixelType, 3u>, ContourDeformationType > MeshSource;
	typedef itk::AddImageFilter<ImageType, ImageType, ImageType > Add;


//	typedef itk::VTKPolyDataReader< ContourDeformationType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDeformationType >     WriterType;
	typedef itk::ImageFileWriter<ImageType>                      ImageWriter;

	typedef itk::EllipseSpatialObject< Dimension >   EllipseType;
	typedef itk::SpatialObjectToImageFilter< EllipseType, ImageType >   SpatialObjectToImageFilterType;

	typedef itk::NormalQuadEdgeMeshFilter
		    < ContourDeformationType, ContourDeformationType > NormalFilterType;
		typedef typename NormalFilterType::Pointer             NormalFilterPointer;


	SpatialObjectToImageFilterType::Pointer imageFilter1 = SpatialObjectToImageFilterType::New();
	ImageType::SizeType size;
	size.Fill( 100 );
	imageFilter1->SetSize( size );
	ImageType::SpacingType spacing; spacing.Fill(1.0);
	imageFilter1->SetSpacing( spacing );

	EllipseType::Pointer ellipse1 = EllipseType::New();
	EllipseType::ArrayType radius; radius.Fill( 20 );
	radius[0] = 10;
	radius[1] = 30;
	ellipse1->SetRadius( radius );

	typedef EllipseType::TransformType TransformType;
	TransformType::Pointer tf = TransformType::New();
	tf->SetIdentity();
	TransformType::OutputVectorType traslation;
	TransformType::CenterType center;

	for(size_t i = 0; i<Dimension; i++) traslation[i] = size[i] * 0.5;
	tf->Translate( traslation, false );
	ellipse1->SetObjectToParentTransform( tf );
	ellipse1->SetDefaultOutsideValue( 0.0 );
	ellipse1->SetDefaultInsideValue( 255.0 );

	ellipse1->Update();

	float inside=127;
	imageFilter1->SetInput( ellipse1 );;
	imageFilter1->SetUseObjectValue( false );
	imageFilter1->SetOutsideValue( 0.0 );
	imageFilter1->SetInsideValue( inside );
	imageFilter1->Update();

	EllipseType::Pointer ellipse2 = EllipseType::New();
	EllipseType::ArrayType radius2;
	for(size_t i = 0; i<Dimension; i++) radius2[i] = radius[i] +10;
    ellipse2->SetRadius( radius2 );
    ellipse2->SetObjectToParentTransform( tf );
    ellipse2->SetDefaultInsideValue(255.0);
    ellipse2->SetDefaultOutsideValue(0.0);
    ellipse2->Update();

	SpatialObjectToImageFilterType::Pointer imageFilter2 = SpatialObjectToImageFilterType::New();
	imageFilter2->SetSize( size );
	imageFilter2->SetSpacing( spacing );
	imageFilter2->SetInput( ellipse2 );
	imageFilter2->SetOutsideValue(0.0);
	imageFilter2->SetInsideValue(128);
	imageFilter2->SetUseObjectValue(false);
	imageFilter2->Update();

	Add::Pointer add = Add::New();
	add->SetInput1( imageFilter1->GetOutput() );
	add->SetInput2( imageFilter2->GetOutput() );
	add->Update();


	ImageWriter::Pointer w = ImageWriter::New();
	std::stringstream ss; ss << TEST_DATA_DIR << "ellipse" << Dimension << "D";
	w->SetFileName( ss.str() + ".nii.gz" );
	w->SetInput(add->GetOutput());
	w->Update();

	MeshSource::Pointer p = MeshSource::New();
	p->SetInput( imageFilter1->GetOutput() );
	p->SetObjectValue( inside );
	p->SetRegionOfInterest( imageFilter1->GetOutput()->GetLargestPossibleRegion( ));
	p->Update();
	ContourDeformationPointer result = p->GetOutput();

	// Compute mesh of normals
	NormalFilterPointer normFilter = NormalFilterType::New();
	normFilter->SetInput( result );
	normFilter->Update();
	ContourDeformationPointer normals = normFilter->GetOutput();

	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( result );
	polyDataWriter->SetFileName( ss.str()+".vtk" );
	polyDataWriter->Update();


	typename ContourDeformationType::PointDataContainerPointer points = normals->GetPointData();
	typename ContourDeformationType::PointDataContainerIterator u_it = points->Begin();

	typename ContourDeformationType::PointType position;
	VectorType norm;
	VectorType zero = itk::NumericTraits<VectorType>::Zero;
	while( u_it != points->End() ) {
		position = normals->GetPoint( u_it.Index() );
		norm = (VectorType) u_it.Value();
		for (size_t i =0; i < norm.Size(); i++) norm[i]*=2.0;
		typename ContourDeformationType::PointType newP = position - norm;
		result->SetPoint( u_it.Index(), newP );
		++u_it;
	}

	WriterType::Pointer polyDataWriter2 = WriterType::New();
	polyDataWriter2->SetInput( result );
	polyDataWriter2->SetFileName( ss.str()+"-small.vtk" );
	polyDataWriter2->Update();

}
