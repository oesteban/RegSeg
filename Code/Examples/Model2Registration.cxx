// --------------------------------------------------------------------------
// File:             Model2Registration.cxx
// Date:             14/12/2012
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





#ifndef DATA_DIR
#define DATA_DIR "../Data/ModelGeneration/Model2/"
#endif

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkVTKPolyDataReader.h>
#include <itkVTKPolyDataWriter.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkComposeImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkDisplacementFieldTransform.h>
#include <itkResampleImageFilter.h>
#include <itkOrientImageFilter.h>
#include "MahalanobisLevelSets.h"
#include "GradientDescentLevelSetsOptimizer.h"

using namespace rstk;

int main(int argc, char *argv[]) {
	typedef itk::Image<itk::Vector<float,1u>, 3u>                ChannelType;
	typedef itk::OrientImageFilter< ChannelType, ChannelType >   OrientFilter;
	typedef MahalanobisLevelSets<ChannelType>                    LevelSetsType;
	typedef LevelSetsType::ContourDeformationType                ContourDeformationType;
	typedef ContourDeformationType::Pointer                      ContourDisplacementFieldPointer;
	typedef LevelSetsType::VectorType                            VectorType;
	typedef LevelSetsType::MeanType                              MeanType;
	typedef LevelSetsType::CovarianceType                        CovarianceType;
	typedef LevelSetsType::DeformationFieldType                  DeformationFieldType;

	typedef GradientDescentLevelSetsOptimizer< LevelSetsType >   Optimizer;
	typedef typename Optimizer::Pointer                          OptimizerPointer;

	typedef itk::VTKPolyDataReader< ContourDeformationType >     ReaderType;
	typedef itk::VTKPolyDataWriter< ContourDeformationType >     WriterType;
	typedef itk::ImageFileReader<ChannelType>                      ImageReader;
	typedef itk::ImageFileWriter<ChannelType>                      ImageWriter;
	typedef itk::ImageFileWriter<DeformationFieldType>           DeformationWriter;

	typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;
	typedef itk::VectorResampleImageFilter
			<DeformationFieldType,DeformationFieldType,double>   DisplacementResamplerType;
	typedef itk::BSplineInterpolateImageFunction
			                <DeformationFieldType>               InterpolatorFunction;
	typedef itk::DisplacementFieldTransform<float, 3u>           TransformType;

	typedef itk::ResampleImageFilter<ChannelType,ChannelType,float>    ResamplerType;



	ImageReader::Pointer r = ImageReader::New();
	r->SetFileName( std::string( DATA_DIR ) + "deformed.nii" );
	r->Update();
	ChannelType::Pointer im = r->GetOutput();
	ChannelType::DirectionType dir = im->GetDirection();

	ChannelType::DirectionType ident; ident.SetIdentity();
	OrientFilter::Pointer orient = OrientFilter::New();
	orient->UseImageDirectionOn();
	orient->SetDesiredCoordinateDirection( ident );
	orient->SetInput( im );
	orient->Update();
	im = orient->GetOutput();
	ChannelType::PointType newOrig; newOrig.Fill(0.0);
	im->SetOrigin( newOrig );

	ImageWriter::Pointer ww = ImageWriter::New();
	ww->SetInput( im );
	ww->SetFileName( "refImage.nii" );
	ww->Update();


	VectorType zero = itk::NumericTraits<VectorType>::Zero;

	ReaderType::Pointer polyDataReader1 = ReaderType::New();
	polyDataReader1->SetFileName( std::string( DATA_DIR ) + "gmwm_surface.vtk" );
	polyDataReader1->Update();
	ContourDisplacementFieldPointer initialContour1 = polyDataReader1->GetOutput();

	typename ContourDeformationType::PointsContainerPointer points1 = initialContour1->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it1 = points1->Begin();

	while( u_it1 != points1->End() ) {
		initialContour1->SetPointData( u_it1.Index(),zero);
		++u_it1;
	}

	ReaderType::Pointer polyDataReader2 = ReaderType::New();
	polyDataReader2->SetFileName( std::string( DATA_DIR ) + "pial_surface.vtk" );
	polyDataReader2->Update();
	ContourDisplacementFieldPointer initialContour2 = polyDataReader2->GetOutput();

	typename ContourDeformationType::PointsContainerPointer points2 = initialContour2->GetPoints();
	typename ContourDeformationType::PointsContainerIterator u_it2 = points2->Begin();

	while( u_it2 != points2->End() ) {
		initialContour2->SetPointData( u_it2.Index(),zero);
		++u_it2;
	}

	// Initialize tissue signatures
	MeanType mean0; // This is WM
	mean0[0] = 0.7806;
	CovarianceType cov0;
	cov0(0,0) =  0.09;
	MeanType mean1; // This is GM
	mean1[0] = 0.3452;
	CovarianceType cov1;
	cov1(0,0) =  0.09;
	MeanType mean3; // This is CSF
	mean3[0] = 0.2363;
	CovarianceType cov3;
	cov3(0,0) =  0.09;

	typename LevelSetsType::ParametersType params1;
	params1.mean[0] = mean0;
	params1.mean[1] = mean1;
	params1.iCovariance[0] = cov0;
	params1.iCovariance[1] = cov1;

	typename LevelSetsType::ParametersType params2;
	params2.mean[0] = mean1;
	params2.mean[1] = mean3;
	params2.iCovariance[0] = cov1;
	params2.iCovariance[1] = cov3;

	// Initialize deformation field
	DeformationFieldType::Pointer df = DeformationFieldType::New();
	DeformationFieldType::SizeType imSize = im->GetLargestPossibleRegion().GetSize();
	DeformationFieldType::SizeType size; size.Fill(16);
	DeformationFieldType::SpacingType spacing = im->GetSpacing();
	for( size_t i = 0; i<3; i++) spacing[i] = 1.0*imSize[i]/size[i];
	df->SetRegions( size );
	df->SetSpacing( spacing );
	df->SetDirection( im->GetDirection() );
	df->Allocate();
	df->FillBuffer( itk::NumericTraits<DeformationFieldType::PixelType>::Zero );

	// Initialize LevelSet function
	LevelSetsType::Pointer ls = LevelSetsType::New();
	ls->SetReferenceImage( im );
	ls->AddShapePrior( initialContour1, params1 );
	ls->AddShapePrior( initialContour2, params2 );

	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetLevelSetsFunction( ls );
	opt->SetDeformationField( df );
	opt->SetNumberOfIterations(500);
	opt->SetAlpha( 50 );
	opt->SetBeta( 100 );
	opt->SetStepSize( 0.001 );

	// Start
	opt->Start();

	// Write final result out
	WriterType::Pointer polyDataWriter = WriterType::New();
	polyDataWriter->SetInput( ls->GetCurrentContourPosition()[0] );
	polyDataWriter->SetFileName( "Model2-wm.vtk" );
	polyDataWriter->Update();

	// Write final result out
	WriterType::Pointer polyDataWriter2 = WriterType::New();
	polyDataWriter2->SetInput( ls->GetCurrentContourPosition()[1] );
	polyDataWriter2->SetFileName( "Model2-pial.vtk" );
	polyDataWriter2->Update();

	DeformationWriter::Pointer w = DeformationWriter::New();
	w->SetInput( opt->GetDeformationField() );
	w->SetFileName( "Model2_field.nii.gz" );
	w->Update();

	DisplacementResamplerType::Pointer p = DisplacementResamplerType::New();
	p->SetInput( opt->GetDeformationField() );
	p->SetOutputOrigin( im->GetOrigin() );
	p->SetOutputSpacing( im->GetSpacing() );
	p->SetOutputDirection( im->GetDirection() );
	p->SetSize( im->GetLargestPossibleRegion().GetSize() );
	//p->SetInterpolator( InterpolatorFunction::New() );
	p->Update();
	DeformationWriter::Pointer w2 = DeformationWriter::New();
	w2->SetInput( p->GetOutput() );
	w2->SetFileName( "Model2_fieldHD.nii.gz" );
	w2->Update();
/*
	DeformationFieldType::Pointer dfield = DeformationFieldType::New();
	TransformType::Pointer tf = TransformType::New();
	tf->SetDisplacementField( dfield );

	ResamplerType::Pointer res = ResamplerType::New();
	res->SetInput( im );
	res->SetReferenceImage( im );
	res->SetTransform( tf );
	res->UseReferenceImageOn();
	res->Update();

	ImageWriter::Pointer w3 = ImageWriter::New();
	w3->SetInput( res->GetOutput() );
	w3->SetFileName( "FAunwarped.nii.gz" );
	w3->Update();*/
}


