// --------------------------------------------------------------------------------------
// File:          SegmentationOptimizer.hxx
// Date:          Feb 10, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// ACWE-Reg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWE-Reg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWE-Reg.  If not, see <http://www.gnu.org/licenses/>.
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
