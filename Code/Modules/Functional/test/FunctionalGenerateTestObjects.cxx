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
#include "MahalanobisFunctional.h"

using namespace rstk;


const unsigned int Dimension = 3;

int main(int argc, char *argv[]) {
	typedef float                                PixelType;
	typedef itk::Vector<float, 2u>               VectorPixelType;
	typedef itk::Image<PixelType, Dimension>            ImageType;
	typedef MahalanobisFunctional<itk::Image<VectorPixelType, Dimension> >      FunctionalType;

	typedef typename FunctionalType::ContourType     ContourType;
	typedef typename ContourType::Pointer           ContourDeformationPointer;
	typedef typename FunctionalType::PointType                 PointType;

	typedef itk::BinaryMask3DMeshSource< itk::Image<PixelType, 3u>, ContourType > MeshSource;
	typedef itk::AddImageFilter<ImageType, ImageType, ImageType > Add;


//	typedef itk::VTKPolyDataReader< ContourType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourType >     WriterType;
	typedef itk::ImageFileWriter<ImageType>                      ImageWriter;

	typedef itk::EllipseSpatialObject< Dimension >   EllipseType;
	typedef itk::SpatialObjectToImageFilter< EllipseType, ImageType >   SpatialObjectToImageFilterType;

	typedef itk::NormalQuadEdgeMeshFilter
		    < ContourType, ContourType > NormalFilterType;
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


	typename ContourType::PointDataContainerPointer points = normals->GetPointData();
	typename ContourType::PointDataContainerIterator u_it = points->Begin();

	typename ContourType::PointType position;
	PointType norm;
	PointType zero = itk::NumericTraits<PointType>::Zero;
	while( u_it != points->End() ) {
		position = normals->GetPoint( u_it.Index() );
		norm = (PointType) u_it.Value();
		for (size_t i =0; i < norm.Size(); i++) norm[i]*=2.0;
		typename ContourType::PointType newP = position - norm;
		result->SetPoint( u_it.Index(), newP );
		++u_it;
	}

	WriterType::Pointer polyDataWriter2 = WriterType::New();
	polyDataWriter2->SetInput( result );
	polyDataWriter2->SetFileName( ss.str()+"-small.vtk" );
	polyDataWriter2->Update();

}
