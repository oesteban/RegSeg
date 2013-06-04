// --------------------------------------------------------------------------------------
// File:          LevelSetsBaseTest.cxx
// Date:          Jun 4, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "gtest/gtest.h"

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "/home/oesteban/workspace/ACWE-Reg/Data/Ellipse/"
// data source http://code.google.com/p/v3dsolver/source/browse/data/
#endif

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkQuadEdgeMesh.h>
#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>
#include <itkVectorImageToImageAdaptor.h>

#include "MahalanobisLevelSets.h"

using namespace rstk;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

typedef itk::Vector<float, 1u>               VectorPixelType;
typedef itk::Image<VectorPixelType, 3u>      ImageType;
typedef MahalanobisLevelSets<ImageType>      LevelSetsType;

typedef LevelSetsType::ContourDeformationType     ContourDeformationType;
typedef ContourDeformationType::Pointer           ContourDisplacementFieldPointer;
typedef LevelSetsType::VectorType                 VectorType;
typedef LevelSetsType::MeanType                   MeanType;
typedef LevelSetsType::CovarianceType             CovarianceType;
typedef LevelSetsType::DeformationFieldType       DeformationFieldType;

typedef itk::MeshFileReader< ContourDeformationType >        ReaderType;
typedef itk::MeshFileWriter< ContourDeformationType >        WriterType;
typedef itk::ImageFileReader<ImageType>                      ImageReader;
typedef itk::ImageFileWriter<ImageType>                      ImageWriter;
typedef rstk::DisplacementFieldFileWriter<DeformationFieldType> Writer;

typedef itk::VectorImageToImageAdaptor<double,3u>            VectorToImage;



class MahalanobisLevelSetsTest : public ::testing::Test {
public:
	virtual void SetUp() {
		ImageReader::Pointer r = ImageReader::New();
		std::string fname1 = std::string( TEST_DATA_DIR ) + "ellipse3D.nii.gz";
		r->SetFileName( fname1 );
		EXPECT_NO_THROW( r->Update() ) << "Failed to load '" << fname1 << "'";

		reference = r->GetOutput();

		ReaderType::Pointer polyDataReader = ReaderType::New();
		std::string fname2 = std::string( TEST_DATA_DIR ) + "ellipse3D.vtk";
		polyDataReader->SetFileName( fname2 );
		EXPECT_NO_THROW( polyDataReader->Update() )  << "Failed to load '" << fname2 << "'";
		prior = polyDataReader->GetOutput();


		ls = LevelSetsType::New();
		ls->SetReferenceImage( reference );
		ls->AddShapePrior( prior );
		ls->Initialize();
	}

	ImageType::Pointer reference;
	ContourDisplacementFieldPointer prior;
	LevelSetsType::Pointer ls;
};


TEST_F( MahalanobisLevelSetsTest, DeformationParametersGridExtent ) {
	ASSERT_TRUE( true );
}
