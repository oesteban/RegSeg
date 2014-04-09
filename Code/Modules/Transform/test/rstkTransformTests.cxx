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
#include <itkResampleImageFilter.h>
#include <itkImageAlgorithm.h>
#include <itkImageFileReader.h>
#include <itkBSplineInterpolateImageFunction.h>
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


typedef float ScalarType;
typedef itk::ContinuousIndex< ScalarType, 3> CIndex;

typedef itk::Point<ScalarType, 3> PointType;
typedef itk::Vector<ScalarType, 3 > VectorType;
typedef itk::Image<ScalarType, 3> ComponentType;
typedef itk::Image< VectorType, 3 > FieldType;
typedef itk::ImageFileReader< ComponentType > FieldReader;
typedef itk::ComposeImageFilter< ComponentType, FieldType > ComposeFilter;
typedef itk::ResampleImageFilter< ComponentType, ComponentType > ExpandFilter;
typedef itk::BSplineInterpolateImageFunction< ComponentType, double, double > BSplineInterpolator;
typedef rstk::DisplacementFieldComponentsFileWriter<FieldType> Writer;

typedef BSplineSparseMatrixTransform<ScalarType, 3, 3> Transform;
typedef Transform::Pointer                             TPointer;
typedef typename Transform::CoefficientsImageType      CoefficientsType;

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
	FieldType::Pointer m_field, m_orig_field, m_hr_field;
	size_t m_N;
	size_t m_K;

	void InitHRField( float factor = 2.0 ) {
		FieldType::DirectionType newDir = m_orig_field->GetDirection();

		FieldType::SizeType     newSize = m_orig_field->GetLargestPossibleRegion().GetSize();
		for ( size_t i = 0; i < 3; i++ ) newSize[i]= floor( factor * newSize[i] ) + 1;

		CIndex start; start.Fill( -0.5 );
		CIndex end;
		for ( size_t i = 0; i < 3; i++ )
			end[i]= m_orig_field->GetLargestPossibleRegion().GetSize()[i] - 0.5;

		PointType domS;
		m_orig_field->TransformContinuousIndexToPhysicalPoint( start, domS );
		PointType domE;
		m_orig_field->TransformContinuousIndexToPhysicalPoint(   end, domE );

		typedef itk::Matrix< typename VectorType::ComponentType, 3, 3 > VectorValuedMatrix;
		VectorValuedMatrix m; m.Fill(0.0);

		for( size_t i = 0; i< 3; i++ )
			for (size_t j=0; j<3; j++ )
				m(i,j) = static_cast<typename VectorType::ComponentType>( newDir(i,j) );

		FieldType::SpacingType newSpacing;
		VectorType extent = domE - domS;
		VectorType oldExt = m * extent;
		VectorType hSpacing;

		for ( size_t i = 0; i < 3; i++ ) {
			newSpacing[i] = fabs( oldExt[i] )/(1.0*newSize[i]);
		}

		PointType newOrigin = domS + newDir * newSpacing * 0.5;


		ComposeFilter::Pointer c = ComposeFilter::New();
		for (size_t i = 0; i<3; i++) {
			std::stringstream ss;
			ss << TEST_DATA_DIR << "field_" << i << "_lr.nii.gz";
			FieldReader::Pointer r = FieldReader::New();
			r->SetFileName( ss.str().c_str() );
			r->Update();

			ExpandFilter::Pointer e = ExpandFilter::New();
			e->SetInput( r->GetOutput() );
			e->SetInterpolator( BSplineInterpolator::New() );
			e->SetOutputDirection( newDir );
			e->SetOutputOrigin( newOrigin );
			e->SetSize( newSize );
			e->SetOutputSpacing( newSpacing );
			e->Update();
			c->SetInput( i, e->GetOutput() );

		}
		c->Update();
		m_hr_field = c->GetOutput();

		Writer::Pointer w = Writer::New();
		w->SetInput( m_hr_field );
		std::stringstream ss;
		ss << "hr_field_" << factor;
		w->SetFileName( ss.str().c_str() );
		w->Update();
	}

};


TEST_F( TransformTests, MatricesTest ) {
	m_transform->ComputeCoefficients();
	m_transform->UpdateField();
	m_transform->SetOutputReference( m_field );
	m_transform->Interpolate();
	ASSERT_TRUE( m_transform->GetPhi() == m_transform->GetS() );
}

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
	w->SetFileName( "orig_field_resampled");
	w->Update();

	const VectorType* rbuf = m_orig_field->GetBufferPointer();
	const VectorType* tbuf = m_field->GetBufferPointer();

	VectorType v1, v2;
	double error = 0.0;
	for( size_t i = 0; i< this->m_K; i++ ) {
		v1 = *( rbuf + i );
		v2 = *( tbuf + i );

		error+= ( v2 - v1 ).GetNorm();
	}
	error = error * (1.0/ m_transform->GetNumberOfParameters() );

	ASSERT_NEAR( 0.0, error, 1.0e-5 );
}

TEST_F( TransformTests, InterpolateOneSample1 ) {
	PointType p;
	FieldType::IndexType idx;
	FieldType::SizeType s = m_field->GetLargestPossibleRegion().GetSize();

	for ( size_t i = 0; i<3; i++ ) {
		idx[i] = floor( (s[i]-1)*0.5 );
	}

	m_field->TransformIndexToPhysicalPoint( idx, p );
	VectorType v1 = m_field->GetPixel( idx );

	m_transform->ComputeCoefficients();
	m_transform->UpdateField();
	m_transform->AddOffGridPos( p );

	m_transform->Interpolate();

	Transform::WeightsMatrix m = m_transform->GetPhi();
	VectorType v2 = m_transform->GetOffGridValue( 0 );

	ASSERT_NEAR( 0, (v1-v2).GetNorm(), 1.0e-5 );

}

TEST_F( TransformTests, InterpolateOneSample2 ) {
	m_transform->ComputeCoefficients();
	m_transform->UpdateField();
	this->InitHRField( 2.0 );

	PointType p;
	FieldType::IndexType idx;
	FieldType::SizeType s = m_hr_field->GetLargestPossibleRegion().GetSize();

	for ( size_t i = 0; i<3; i++ ) {
		idx[i] = floor( (s[i]-1)*0.5 );
	}
	m_hr_field->TransformIndexToPhysicalPoint( idx, p );
	m_transform->AddOffGridPos( p );
	m_transform->Interpolate();

	Transform::WeightsMatrix m = m_transform->GetPhi();

	VectorType v1 = m_hr_field->GetPixel( idx );
	VectorType v2 = m_transform->GetOffGridValue( 0 );
	ASSERT_NEAR( 0, (v1-v2).GetNorm(), 1.0e-5 );

}

TEST_F( TransformTests, InterpolateOneSample3 ) {
	m_transform->ComputeCoefficients();
	m_transform->UpdateField();
	this->InitHRField( 2.3 );

	PointType p;
	FieldType::IndexType idx;
	FieldType::SizeType s = m_hr_field->GetLargestPossibleRegion().GetSize();

	for ( size_t i = 0; i<3; i++ ) {
		idx[i] = floor( (s[i]-1)*0.5 );
	}
	m_hr_field->TransformIndexToPhysicalPoint( idx, p );
	m_transform->AddOffGridPos( p );
	m_transform->Interpolate();

	VectorType v1 = m_hr_field->GetPixel( idx );
	VectorType v2 = m_transform->GetOffGridValue( 0 );
	ASSERT_NEAR( 0, (v1-v2).GetNorm(), 1.0e-5 );

}

TEST_F( TransformTests, InterpolateAllSamples1 ) {
	m_transform->ComputeCoefficients();
	m_transform->UpdateField();

	this->InitHRField( 2.3 );
	size_t nSamples =  m_hr_field->GetLargestPossibleRegion().GetNumberOfPixels();

	PointType p;
	FieldType::IndexType idx;

	for ( size_t i = 0; i<nSamples; i++ ) {
		idx = m_hr_field->ComputeIndex( i );
		m_hr_field->TransformIndexToPhysicalPoint( idx, p );
		m_transform->AddOffGridPos( p );
	}

	m_transform->Interpolate();

	Transform::WeightsMatrix m = m_transform->GetPhi();
	Transform::SparseVectorType row;

	Transform::CoefficientsImageArray coeff = m_transform->GetCoefficientsImages();
	size_t c;

	VectorType v1;
	VectorType v2;
	double wi, coeffk;
	double error = 0.0;
	for( size_t r = 0; r<m.rows(); r++ ) {
		row = m.get_row( r );
		v2 = m_transform->GetOffGridValue( r );
		v1.Fill( 0.0 );

		for( size_t i = 0; i<row.size(); i++ ) {
			c = row[i].first;
			wi = row[i].second;

			for( size_t j = 0; j<3; j++ ) {
				coeffk = coeff[j]->GetPixel( coeff[j]->ComputeIndex(c) );
				v1[j]+= coeffk * wi;
			}
		}
		error+= (v1-v2).GetNorm();
	}

	error/=nSamples;

	ASSERT_NEAR( 0, error, 1.0e-5 );

}

TEST_F( TransformTests, InterpolateAllSamples2 ) {
	m_transform->ComputeCoefficients();
	m_transform->UpdateField();

	this->InitHRField( 2.0 );
	size_t nSamples =  m_hr_field->GetLargestPossibleRegion().GetNumberOfPixels();

	PointType p;
	FieldType::IndexType idx;

	for ( size_t i = 0; i<nSamples; i++ ) {
		idx = m_hr_field->ComputeIndex( i );
		m_hr_field->TransformIndexToPhysicalPoint( idx, p );
		m_transform->AddOffGridPos( p );
	}

	m_transform->Interpolate();


	VectorType v1, v2;
	double error = 0.0;
	for ( size_t i = 0; i<nSamples; i++ ) {
		idx = m_hr_field->ComputeIndex( i );
		v1 = m_hr_field->GetPixel( idx );
		v2 = m_transform->GetOffGridValue( i );
		error+= (v1-v2).GetNorm();
	}
	error/=nSamples;
	ASSERT_NEAR( 0, error, 1.0e-5 );

}


TEST_F( TransformTests, CompareBSplineInterpolation ) {
	this->InitHRField( 2.0 );
	m_transform->ComputeCoefficients();
	m_transform->SetOutputReference( m_hr_field );
	m_transform->Interpolate();

	Writer::Pointer w = Writer::New();
	w->SetInput( m_transform->GetOutputField() );
	w->SetFileName( "interpolated_field");
	w->Update();

	const VectorType* rbuf = m_hr_field->GetBufferPointer();
	const VectorType* tbuf = m_transform->GetOutputField()->GetBufferPointer();

	VectorType v1, v2;
	double error = 0.0;
	for( size_t i = 0; i< m_transform->GetNumberOfSamples(); i++ ) {
		v1 = *( rbuf + i );
		v2 = *( tbuf + i );
		error+= ( v2 - v1 ).GetNorm();
	}
	error*= (1.0/ m_transform->GetNumberOfSamples() );
	ASSERT_NEAR( 0.0, error, 1.0e-1 );

}


TEST_F( TransformTests, RandomSampleTest ) {
	this->InitHRField( 3.27 );

	m_transform->ComputeCoefficients();
	m_transform->UpdateField();

	m_transform->SetOutputReference( m_hr_field );
	m_transform->Interpolate();

	Writer::Pointer w = Writer::New();
	w->SetInput( m_transform->GetOutputField() );
	w->SetFileName( "interpolated_field_3.7");
	w->Update();


	TPointer tfm = Transform::New();
	tfm->CopyGridInformation( m_transform->GetCoefficientsImages()[0] );
	tfm->SetCoefficientsImages( m_transform->GetCoefficientsImages() );

	int rindex = 0;
	std::vector< int > subsample;
	Transform::PointType p;
	for ( size_t i = 0; i<200; i++) {
		rindex = rand() % (m_transform->GetNumberOfSamples() + 1);
		subsample.push_back( rindex );
		m_hr_field->TransformIndexToPhysicalPoint( m_hr_field->ComputeIndex(rindex) ,p );
		tfm->AddOffGridPos( p );
	}

	tfm->Interpolate();

	FieldType::ConstPointer field = m_transform->GetOutputField();

	double error = 0.0;
	VectorType v1, v2;
	for ( size_t i = 0; i<subsample.size(); i++ ) {
		v1 = tfm->GetOffGridValue( i );
		v2 = field->GetPixel( m_hr_field->ComputeIndex( subsample[i] ) );

		//EXPECT_NEAR( v2.GetNorm(), v1.GetNorm(), 1.0e-3 );

		error+= (v1-v2).GetNorm();
	}
	ASSERT_NEAR( 0.0, error, 1.0e-5 );
}
} // namespace rstk
