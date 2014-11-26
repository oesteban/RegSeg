// --------------------------------------------------------------------------------------
// File:          CompositeMatrixTransform.hxx
// Date:          Oct 30, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

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
::InitializeField() {
	if( this->m_DisplacementField.IsNull() ) {
		OutputVectorType zerov; zerov.Fill(0.0);
		this->m_DisplacementField = DisplacementFieldType::New();
		this->m_DisplacementField->SetRegions(this->m_ReferenceSize);
		this->m_DisplacementField->SetSpacing(this->m_ReferenceSpacing);
		this->m_DisplacementField->SetOrigin(this->m_ReferenceOrigin);
		this->m_DisplacementField->SetDirection(this->m_DomainDirection);
		this->m_DisplacementField->Allocate();
		this->m_DisplacementField->FillBuffer(zerov);
	}
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::Compute() {
	if ((this->m_NumberOfTransforms == 0) || (this->m_Components.size()!=this->m_NumberOfTransforms)) {
		itkExceptionMacro(<< "number of transforms is zero or it does not match the number of stored coefficients sets.");
	}

	if( this->m_InterpolationMode == Superclass::UNKNOWN ) {
		itkExceptionMacro(<< "output mode has not been initialized");
	}
	this->InitializeField();

	DisplacementType* dfbuffer = this->m_DisplacementField->GetBufferPointer();
	size_t nPix = this->m_DisplacementField->GetLargestPossibleRegion().GetNumberOfPixels();

	DisplacementType v;
	for( size_t c = 0; c < this->m_NumberOfTransforms; c++) {
		//itkDynamicCastInDebugMode< const CoefficientsImageType* >
		this->m_Components[c]->Interpolate();
		FieldConstPointer f = this->m_Components[c]->GetField();

		const DisplacementType* cbuffer = f->GetBufferPointer();
		DisplacementType vc;
		for( size_t i = 0; i < nPix; i++) {
			vc = *(cbuffer + i);
			v = *(dfbuffer + i);
			if( vc.GetNorm() > 1.e-5)
				*(dfbuffer + i)+= vc;
		}
	}
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::PushBackCoefficients(const CoefficientsImageArray & images) {
	BSplineComponentPointer tf = BSplineComponentType::New();

	if(this->m_InterpolationMode == Superclass::POINTS_MODE) {
		tf->SetOutputPoints( this->GetPointLocations() );
	} else if (this->m_InterpolationMode == Superclass::GRID_MODE) {
		tf->SetOutputReference(this->m_DisplacementField);
	}

	tf->SetCoefficientsImages(images);

	this->PushBackTransform(tf);
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::PushBackCoefficients(const FieldType* field) {
	BSplineComponentPointer tf = BSplineComponentType::New();
	if(this->m_InterpolationMode == Superclass::POINTS_MODE) {
		tf->SetOutputPoints( this->GetPointLocations() );
	} else if (this->m_InterpolationMode == Superclass::GRID_MODE) {
		tf->SetOutputReference(this->m_DisplacementField);
	}

	tf->SetCoefficientsVectorImage(field);
	this->PushBackTransform(tf);
}

} // namespace rstk
