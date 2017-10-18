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

#ifndef SEGMENTATIONOPTIMIZER_HXX_
#define SEGMENTATIONOPTIMIZER_HXX_


#include "SegmentationOptimizer.h"

#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <itkComplexToRealImageFilter.h>
#include <itkImageAlgorithm.h>

using namespace std;

#include "DisplacementFieldFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"

namespace rstk {

/**
 * Default constructor
 */
template< typename TFunctional >
SegmentationOptimizer<TFunctional>::SegmentationOptimizer(): Superclass() {}

template< typename TFunctional >
void SegmentationOptimizer<TFunctional>
::PrintSelf(std::ostream &os, itk::Indent indent) const {
	Superclass::PrintSelf(os,indent);
	os << indent << "Learning rate:" << this->m_LearningRate << std::endl;
	os << indent << "Number of iterations: " << this->m_NumberOfIterations << std::endl;
	os << indent << "Current iteration: " << this->m_CurrentIteration << std::endl;
	os << indent << "Stop condition:" << this->m_StopCondition << std::endl;
	os << indent << "Stop condition description: " << this->m_StopConditionDescription.str() << std::endl;
}


template< typename TFunctional >
void SegmentationOptimizer<TFunctional>::InitializeParameters() {

}

template< typename TFunctional >
void SegmentationOptimizer<TFunctional>::ComputeDerivative() {
	// Multiply phi and copy reshaped on this->m_Derivative
	this->m_Gradient = this->m_Functional->ComputeDerivative();

	if ( this->m_Displacement.rows() == 0 || this->m_Displacement.cols() == 0) {
		this->m_Displacement = WeightsMatrix( this->m_Gradient.rows(), this->m_Gradient.cols() );
	}
}

template< typename TFunctional >
void SegmentationOptimizer<TFunctional>::Iterate() {
	SampleType sample;
	size_t nrows = this->m_Gradient.rows();
	PointValueType g;
	VectorType ci;
	PointValueType gradSum = 0.0;

	for ( size_t i = 0; i < nrows; i++ ) {
		g = 0.0;
		ci.Fill( 0.0 );

		if( !this->m_Gradient.empty_row(i) ) {
			for( size_t d = 0; d < Dimension; d++) {
				ci[d] = this->m_Gradient( i, d );
			}
			g = ci.GetNorm();
			gradSum+=g;
		}

		sample.push_back( GradientSample( g, ci, i ) );
	}

	std::sort(sample.begin(), sample.end());
	size_t sSize = sample.size();
	size_t q1 = floor( (sSize-1)* 0.1);
	size_t q2 = round( (sSize-1)*0.50 );
	size_t q3 = ceil ( (sSize-1)* (1.0 - 0.65 ) );

#ifndef NDEBUG
	std::cout << "\tavg=" << (gradSum/sSize) << ", max=" << sample[sSize-1].grad << ", min=" << sample[0].grad << ", q1=" << sample[q1].grad << ", q2=" << sample[q3].grad << ", med=" << sample[q2].grad << "." << std::endl;
#endif

	VectorType vi;
	PointValueType old_grad;
	vnl_random rnd = vnl_random();
	for( size_t i = q3; i<sample.size(); i++ ){
		old_grad = sample[i].grad;
		ci.Fill(0.0);
		vi.Fill(0.0);
		if (old_grad>0.0) {
			ci = (sample[i].normal) / old_grad;
			g = sample[rnd.lrand32(0,q3-1)].grad;
			sample[i].grad = g;
			vi = ci*g;
			sample[i].normal = vi;
		}

		for( size_t j = 0; j< Dimension; j++ ) {
			this->m_Gradient.put( sample[i].gid, j, vi[j] );
		}
	}

#ifndef NDEBUG
	std::sort(sample.begin(), sample.end() );
	std::cout << "\tavg=" << (gradSum/sSize) << ", max=" << sample[sSize-1].grad << ", min=" << sample[0].grad << ", q1=" << sample[q1].grad << ", q2=" << sample[q3].grad << ", med=" << sample[q2].grad << "." << std::endl;
#endif

}

template< typename TFunctional >
void SegmentationOptimizer<TFunctional>::PostIteration() {
	PointValueType g;
	VectorType ci;
	PointValueType gradSum = 0.0;
	size_t nrows = this->m_Gradient.rows();
	for ( size_t i = 0; i < nrows; i++ ) {
		if( !this->m_Gradient.empty_row(i) ) {
			ci.Fill(0.0);
			for( size_t d = 0; d < Dimension; d++) {
				ci[d] = this->m_Gradient( i, d );
				this->m_Displacement(i,d) += ci[d];
			}
			g = ci.GetNorm();
			gradSum+=g;
		}
	}

	this->m_Functional->SetCurrentDisplacements( this->m_Displacement );
}


} // end namespace rstk




#endif /* SEGMENTATIONOPTIMIZER_HXX_ */
