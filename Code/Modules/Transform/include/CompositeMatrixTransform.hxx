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

#include "CompositeMatrixTransform.h"
#include "rstkCoefficientsWriter.h"

namespace rstk {

template< class TScalar, unsigned int NDimensions >
CompositeMatrixTransform<TScalar,NDimensions>
::CompositeMatrixTransform():
Superclass(),
m_NumberOfTransforms(0) {}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::PrintSelf(std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
	os << indent << indent << "Number of transforms: "<< this->m_NumberOfTransforms << std::endl;
}


template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::Interpolate() {
	if ((this->m_NumberOfTransforms == 0) || (this->m_Components.size()!=this->m_NumberOfTransforms)) {
		itkExceptionMacro(<< "number of transforms is zero or it does not match the number of stored coefficients sets.");
	}

	switch(this->m_InterpolationMode) {
	case Superclass::GRID_MODE:
		this->ComputeGrid();
		break;
	case Superclass::POINTS_MODE:
		this->ComputePoints();
		break;
	default:
		itkExceptionMacro(<< "output mode has not been initialized");
		break;
	}

}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::ConsistencyCheck() {
	for( size_t c = 0; c < this->m_NumberOfTransforms; c++) {
		if( this->m_Components[c]->GetInterpolationMode() != this->m_InterpolationMode) {
			TransformComponentPointer tf = this->m_Components[c]->Clone();
			this->m_Components[c] = tf;
		}
	}

}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::ComputeGrid() {
	DisplacementType* dfbuffer = this->m_DisplacementField->GetBufferPointer();
	size_t nPix = this->m_DisplacementField->GetLargestPossibleRegion().GetNumberOfPixels();

	DisplacementType vc;
	for( size_t c = 0; c < this->m_NumberOfTransforms; c++) {
		this->m_Components[c]->InterpolateField();
		FieldConstPointer f = this->m_Components[c]->GetDisplacementField();

		const DisplacementType* cbuffer = f->GetBufferPointer();
		for( size_t i = 0; i < nPix; i++) {
			vc = *(cbuffer + i);
			if( vc.GetNorm() > 1.e-5)
				*(dfbuffer + i)+= vc;
		}
	}
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::ComputePoints() {
	size_t nps = this->m_PointLocations.size();

	for( size_t d = 0; d<Dimension; d++) {
		this->m_PointValues[d] = DimensionVector();
		this->m_PointValues[d].set_size( this->m_NumberOfPoints );
		this->m_PointValues[d].fill(0.0);
	}

	VectorType vc;
	for( size_t c = 0; c < this->m_NumberOfTransforms; c++) {
		this->m_Components[c]->SetOutputPoints( this->GetPointLocations() );
		this->m_Components[c]->InterpolatePoints();

		DimensionParameters disp = this->m_Components[c]->GetPointValues();

		for(size_t d = 0; d < Dimension; d++) {
			this->m_PointValues[d]+= disp[d];
		}
	}
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::PushBackCoefficients(const CoefficientsImageArray & images) {
	BSplineComponentPointer tf = BSplineComponentType::New();
	tf->SetCoefficientsImages(images);
	this->PushBackTransform(tf);
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::PushBackCoefficients(const FieldType* field) {
	BSplineComponentPointer tf = BSplineComponentType::New();
	tf->SetCoefficientsVectorImage(field);
	this->PushBackTransform(tf);
}

} // namespace rstk
