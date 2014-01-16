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
#include <itkVectorImage.h>
#include <itkComposeImageFilter.h>
#include <itkImageAlgorithm.h>
#include <itkImageFileReader.h>
#include "BSplineSparseMatrixTransform.h"
#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "../../Data/Tests/"
#endif

using namespace rstk;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


typedef double ScalarType;
typedef itk::Point<ScalarType, 3> PointType;
typedef itk::Vector<ScalarType, 3 > VectorType;
typedef itk::Image<ScalarType, 3> ComponentType;
typedef itk::Image< VectorType, 3 > FieldType;
typedef itk::ImageFileReader< ComponentType > FieldReader;
typedef itk::ComposeImageFilter< ComponentType, FieldType > ComposeFilter;
typedef BSplineSparseMatrixTransform<ScalarType, 3, 3> Transform;
typedef Transform::Pointer                   TPointer;
typedef rstk::DisplacementFieldComponentsFileWriter<FieldType> Writer;

typedef typename Transform::CoefficientsImageType CoefficientsType;
typedef itk::ImageFileWriter< CoefficientsType >  CoefficientsWriterType;
typedef CoefficientsWriterType::Pointer           CoefficientsWriterPointer;

namespace rstk {

class TransformTests : public ::testing::Test {
public:
	virtual void SetUp() {

		ComposeFilter::Pointer c = ComposeFilter::New();

		for (size_t i = 0; i<3; i++) {
			std::stringstream ss;
			ss << TEST_DATA_DIR << "field_" << i << "_lr.nii.gz";
			FieldReader::Pointer r = FieldReader::New();
			r->SetFileName( ss.str().c_str() );
			r->Update();
			c->SetInput( i, r->GetOutput() );
		}
		c->Update();



		m_K = c->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
		m_orig_field = c->GetOutput();

		m_field = FieldType::New();
		m_field->SetRegions( m_orig_field->GetLargestPossibleRegion() );
		m_field->SetSpacing( m_orig_field->GetSpacing() );
		m_field->SetOrigin( m_orig_field->GetOrigin() );
		m_field->SetDirection( m_orig_field->GetDirection() );
		m_field->Allocate();

		itk::ImageAlgorithm::Copy< FieldType,FieldType >(
				m_orig_field, m_field,
				m_orig_field->GetLargestPossibleRegion(),
				m_field->GetLargestPossibleRegion()
		);


		m_transform = Transform::New();
		m_transform->CopyGridInformation( m_field );
		m_transform->SetField( m_field );
	}

	TPointer m_transform;
	FieldType::Pointer m_field, m_orig_field;
	size_t m_N;
	size_t m_K;

};

TEST_F( TransformTests, SparseMatrixComputeCoeffsTest ) {
	Writer::Pointer w = Writer::New();
	w->SetInput( m_field );
	w->SetFileName( "orig_field");
	w->Update();

	m_transform->ComputeCoefficients();
	for (size_t i = 0; i<3; i++ ) {
		std::stringstream ss;
		ss << "coefficients_" << i << ".nii.gz";
		CoefficientsWriterPointer ww = CoefficientsWriterType::New();
		ww->SetInput( m_transform->GetCoefficientsImages()[i] );
		ww->SetFileName( ss.str().c_str() );
		ww->Update();
	}

	m_transform->UpdateField();

	w->SetInput( m_field );
	w->SetFileName( "updated_field");
	w->Update();

	const VectorType* rbuf = m_orig_field->GetBufferPointer();
	const VectorType* tbuf = m_field->GetBufferPointer();

	VectorType v1, v2;
	double error;
	for( size_t i = 0; i< this->m_K; i++ ) {
		v1 = *( rbuf + i );
		v2 = *( tbuf + i );

		error+= ( v2 - v1 ).GetNorm();
		EXPECT_NEAR( v2.GetNorm(), v1.GetNorm(), 1.0e-3 );
	}
	error = error * (1.0/this->m_K);

	ASSERT_TRUE( error < 1.0e-3 );
}


TEST_F( TransformTests, SparseMatrixInterpolate ) {
	m_transform->ComputeCoefficients();

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

	FieldType::SizeType size;
	FieldType::SpacingType spacing;

	for( size_t d = 0; d<3; d++ ) {
		size[d] = (int) (2.75 * refSize[d]);
		spacing[d] = (refEnd[d] - refOrigin[d])/ ((double) size[d] );
	}

	FieldType::Pointer densefield = FieldType::New();
	densefield->SetOrigin(    refOrigin );
	densefield->SetDirection( m_field->GetDirection() );
	densefield->SetRegions( size );
	densefield->SetSpacing( spacing );
	densefield->Allocate();
	densefield->FillBuffer( zero );

	m_N = densefield->GetLargestPossibleRegion().GetNumberOfPixels();
	m_transform->SetNumberOfSamples(m_N);

	PointType p;
	for( size_t i = 0; i<m_N; i++ ) {
		densefield->TransformIndexToPhysicalPoint( densefield->ComputeIndex(i) ,p );
		m_transform->SetOffGridPos(i,p);
	}

	m_transform->Interpolate();

	VectorType v;
	for( size_t i = 0; i<m_N; i++ ){
		v = m_transform->GetOffGridValue(i);
		densefield->SetPixel( densefield->ComputeIndex(i), v);
	}

	Writer::Pointer w = Writer::New();
	w->SetInput( densefield );
	w->SetFileName( "interpolated_field");
	w->Update();


	ASSERT_TRUE( true );
}

TEST_F( TransformTests, SparseMatrixForwardIDWDumbTest ) {

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
		size[d] = (int) (3.27 * refSize[d]);
		spacing[d] = (refEnd[d] - refOrigin[d])/ ((double) size[d] );
	}


	densefield->SetRegions( size );
	densefield->SetSpacing( spacing );
	densefield->Allocate();
	densefield->FillBuffer( zero );

	m_N = densefield->GetLargestPossibleRegion().GetNumberOfPixels();
	m_transform->SetNumberOfSamples(m_N);

	PointType p;
	for( size_t i = 0; i<m_N; i++ ) {
		densefield->TransformIndexToPhysicalPoint( densefield->ComputeIndex(i) ,p );
		m_transform->SetOffGridPos(i,p);
	}


	m_transform->Interpolate();

	VectorType v;
	for( size_t i = 0; i<m_N; i++ ){
		v = m_transform->GetOffGridValue(i);
		densefield->SetPixel( densefield->ComputeIndex(i), v);
	}

	Writer::Pointer w = Writer::New();
	w->SetInput( densefield );
	w->SetFileName( "dumb_test");
	w->Update();


	m_transform->UpdateField();
	w->SetInput( m_field );
	w->SetFileName( "dumb_test_low" );
	w->Update();

	TPointer tfm = Transform::New();
	tfm->SetNumberOfSamples( 200 );
	tfm->CopyGridInformation( m_transform->GetCoefficientsImages()[0] );
	tfm->SetCoefficientsImages( m_transform->GetCoefficientsImages() );

	int rindex = 0;
	std::vector< int > subsample;
	for ( size_t i = 0; i<200; i++) {
		rindex = rand() % (m_N + 1);
		subsample.push_back( rindex );
		densefield->TransformIndexToPhysicalPoint( densefield->ComputeIndex(rindex) ,p );
		tfm->SetOffGridPos( i, p );
	}

	tfm->Interpolate();


	bool isCorrect = true;
	VectorType v1, v2;
	for ( size_t i = 0; i<subsample.size(); i++ ) {
		v1 = tfm->GetOffGridValue( i );
		v2 = densefield->GetPixel( densefield->ComputeIndex( subsample[i] ) );

		isCorrect = ( (v1-v2).GetNorm() < 1.0e-5 ) &&  isCorrect;
	}


	ASSERT_TRUE( isCorrect );
}

}
