/*
 * rstkTransformTests.cxx
 *
 *  Created on: Jun 10, 2013
 *      Author: oesteban
 */

#include "gtest/gtest.h"

#include <itkPoint.h>
#include <itkVector.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include "SparseMatrixTransform.h"
#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "../Data/Tests/"
#endif

using namespace rstk;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


typedef double ScalarType;
typedef itk::Point<ScalarType, 3> PointType;
typedef itk::Vector<ScalarType, 3 > VectorType;
typedef itk::Image< VectorType, 3 > FieldType;
typedef itk::ImageFileReader< FieldType > FieldReader;
typedef SparseMatrixTransform<ScalarType, 3> Transform;
typedef Transform::Pointer                   TPointer;
typedef rstk::DisplacementFieldComponentsFileWriter<FieldType> Writer;

namespace rstk {

class TransformTests : public ::testing::Test {
public:
	virtual void SetUp() {
		FieldReader::Pointer r = FieldReader::New();
		r->SetFileName( std::string( TEST_DATA_DIR ) + "field.nii.gz" );
		r->Update();

		Writer::Pointer w = Writer::New();
		w->SetInput( r->GetOutput() );
		w->SetFileName( "orig_field.nii.gz");
		w->Update();

		m_field = r->GetOutput();
		m_K = m_field->GetLargestPossibleRegion().GetNumberOfPixels();

		m_Sigma[0] = 8.1;
		m_Sigma[1] = 8.1;
		m_Sigma[2] = 8.1;


		m_transform = Transform::New();
		m_transform->SetK( m_K );
		m_transform->SetSigma( m_Sigma );


		const VectorType* buffer = m_field->GetBufferPointer();
		PointType p;
		for( size_t i = 0; i<m_K; i++ ) {
			m_field->TransformIndexToPhysicalPoint( m_field->ComputeIndex(i) ,p );
			m_transform->SetNode    (i,p);
			m_transform->SetNodeData(i,*( buffer + i));
		}
	}

	TPointer m_transform;
	FieldType::Pointer m_field;
	size_t m_N;
	size_t m_K;
	double m_Sigma[3];

};


TEST_F( TransformTests, SparseMatrixForwardIDWTransformTest ) {

	VectorType zero;
	zero.Fill(0.0);

	FieldType::SizeType refSize = m_field->GetLargestPossibleRegion().GetSize();
	PointType refOrigin = m_field->GetOrigin();
	PointType refEnd;
	FieldType::IndexType lastIdx;
	for( size_t i = 0; i<3; i++ ) {
		lastIdx[i] = refSize[i] - 1;
	}
	m_field->TransformIndexToPhysicalPoint( lastIdx, refEnd );


	FieldType::Pointer densefield = FieldType::New();
	densefield->SetOrigin(    refOrigin );
	densefield->SetDirection( m_field->GetDirection() );


	FieldType::SizeType size;
	FieldType::SpacingType spacing;

	for( size_t d = 0; d<3; d++ ) {
		size[d] = (int) (2.75 * refSize[d]);
		spacing[d] = (refEnd[d] - refOrigin[d])/ ((double) size[d] );
	}


	densefield->SetRegions( size );
	densefield->SetSpacing( spacing );
	densefield->Allocate();
	densefield->FillBuffer( zero );

	m_N = densefield->GetLargestPossibleRegion().GetNumberOfPixels();
	m_transform->SetN(m_N);

	const VectorType* fbuf = m_field->GetBufferPointer();

	PointType p;
	for( size_t i = 0; i<m_N; i++ ) {
		densefield->TransformIndexToPhysicalPoint( densefield->ComputeIndex(i) ,p );
		m_transform->SetPoint    (i,p);
	}


	m_transform->Interpolate();

	VectorType v;
	for( size_t i = 0; i<m_N; i++ ){
		v = m_transform->GetPointData(i);
		densefield->SetPixel( densefield->ComputeIndex(i), v);
	}

	Writer::Pointer w = Writer::New();
	w->SetInput( densefield );
	w->SetFileName( "interpolated_field.nii.gz");
	w->Update();


	ASSERT_TRUE( true );
}


//TEST_F( TransformTests, SparseMatrixBackwardIDWTransformTest ) {
//	VectorType v[3];
//	bool is_equal = true;
//
//	m_transform->ComputeNodes();
//	m_transform->ComputePoints();
//
//	VectorType v2;
//	for( size_t i = 0; i<m_K; i++ ){
//		v2 = m_transform->GetCoefficient(i);
//		if (v2.GetNorm() > 0.0 ) {
//			m_field->SetPixel( m_field->ComputeIndex(i), v2);
//		}
//	}
//
//	Writer::Pointer w = Writer::New();
//	w->SetInput( m_field );
//	w->SetFileName( "coefficients_field.nii.gz");
//	w->Update();
//
//	for( size_t i = 0; i<m_N; i++ ){
//		v[i] = m_transform->GetPointData(i);
//		std::cout << v[i] << " vs. " << this->m_vectors[i] << std::endl;
//	}
//
//	ASSERT_TRUE( true );
//}

}
