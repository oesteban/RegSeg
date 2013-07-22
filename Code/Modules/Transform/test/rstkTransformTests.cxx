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
#include "SparseMatrixTransform.h"
#include "DisplacementFieldFileWriter.h"

using namespace rstk;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


typedef double ScalarType;
typedef itk::Point<ScalarType, 3> PointType;
typedef itk::Vector<ScalarType, 3 > VectorType;
typedef itk::Image< VectorType, 3 > FieldType;
typedef SparseMatrixTransform<ScalarType, 3> Transform;
typedef Transform::Pointer                   TPointer;
typedef rstk::DisplacementFieldFileWriter<FieldType> Writer;

namespace rstk {

class TransformTests : public ::testing::Test {
public:
	virtual void SetUp() {
		m_N = 3;
		m_Sigma[0] = 2.1;
		m_Sigma[1] = 2.1;
		m_Sigma[2] = 2.1;

		m_field = FieldType::New();
		double origin[3] = { 0.0, 0.0, 0.0 };
		FieldType::SizeType size;
		size.Fill( 10 );
		FieldType::DirectionType dir;
		dir.SetIdentity();

		VectorType zero;
		zero.Fill( 0.0 );

		m_field->SetOrigin( origin );
		m_field->SetRegions( size );
		m_field->SetDirection( dir );
		m_field->SetSpacing( 2.0 );

		m_field->Allocate();
		m_field->FillBuffer( zero );

		m_K = m_field->GetLargestPossibleRegion().GetNumberOfPixels();

		m_transform = Transform::New();
		m_transform->SetN( m_N );
		m_transform->SetK( m_K );
		m_transform->SetSigma( m_Sigma );

		m_cp[0][0] = 5.7;
		m_cp[0][1] = 12.2;
		m_cp[0][2] = 3.4;
		m_vectors[0][0] = -10.5;
		m_vectors[0][1] = 7.0;
		m_vectors[0][2] = -3.0;

		m_cp[1][0] = 15.7;
		m_cp[1][1] = 2.2;
		m_cp[1][2] = 3.4;
		m_vectors[1][0] = 15.0;
		m_vectors[1][1] = 7.0;
		m_vectors[1][2] = 3.0;

		m_cp[2][0] = 3.7;
		m_cp[2][1] = 5.2;
		m_cp[2][2] = 13.4;
		m_vectors[2][0] = -12.0;
		m_vectors[2][1] = 27.0;
		m_vectors[2][2] = 13.0;

		for( size_t i = 0; i<m_N; i++ ) {
			m_transform->SetControlPoint    (i,m_cp[i]);
			m_transform->SetControlPointData(i,m_vectors[i]);
		}

		VectorType* fbuf = m_field->GetBufferPointer();

		PointType p;
		FieldType::IndexType idx;
		for( size_t i = 0; i<m_K; i++ ){
			m_field->TransformIndexToPhysicalPoint( m_field->ComputeIndex(i) ,p );
			m_transform->SetGridPoint(i,p);
		}
	}

	TPointer m_transform;
	FieldType::Pointer m_field;
	PointType m_cp[3];
	VectorType m_vectors[3];
	size_t m_N;
	size_t m_K;
	double m_Sigma[3];

};


TEST_F( TransformTests, SparseMatrixForwardIDWTransformTest ) {
	m_transform->ComputeGridPoints();

	VectorType v;
	for( size_t i = 0; i<m_K; i++ ){
		v = m_transform->GetGridPointData(i);
		if (v.GetNorm() > 1.0e-3 ) {
			m_field->SetPixel( m_field->ComputeIndex(i), v);
		}
	}

	Writer::Pointer w = Writer::New();
	w->SetInput( m_field );
	w->SetFileName( "interpolated_field.nii.gz");
	w->Update();


	ASSERT_TRUE( true );
}


TEST_F( TransformTests, SparseMatrixBackwardIDWTransformTest ) {
	VectorType v[3];
	bool is_equal = true;

	m_transform->ComputeGridPoints();

	for( size_t i = 0; i<m_N; i++ ){
		v[i] = m_transform->GetControlPointData(i);
		std::cout << v[i] << " vs. " << this->m_vectors[i] << std::endl;
	}

	for( size_t i = 0; i<m_K; i++ ){
		m_transform->SetGridPointData(i,m_field->GetPixel( m_field->ComputeIndex(i)));
	}
	m_transform->ComputeControlPoints();

	for( size_t i = 0; i<m_N; i++ ){
		v[i] = m_transform->GetControlPointData(i);
		std::cout << v[i] << " vs. " << this->m_vectors[i] << std::endl;
	}

	ASSERT_TRUE( true );
}

}
