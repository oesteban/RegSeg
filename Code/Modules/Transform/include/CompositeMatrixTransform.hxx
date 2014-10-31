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
m_NumberOfTransforms(0) {
	m_ReferenceSize.Fill(0.0);
	m_ReferenceDirection.SetIdentity();
	m_ReferenceOrigin.Fill(0.0);
	m_ReferenceSpacing.Fill(0.0);

}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::PrintSelf(std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);
	os << indent << indent << "Number of transforms: "<< this->m_NumberOfTransforms << std::endl;
	os << indent << indent << "Reference size = " << this->m_ReferenceSize << std::endl;
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::Compute() {
	if ((this->m_NumberOfTransforms == 0) || (this->m_Coefficients.size()!=this->m_NumberOfTransforms)) {
		itkExceptionMacro(<< "number of transforms is zero or it does not match the number of stored coefficients sets.");
	}

	for (size_t d = 0; d < Dimension; d++) {
		if (this->m_ReferenceSize[d] ==0)
			itkExceptionMacro(<< "physical domain is not initialized.");
	}

	DisplacementType zerov; zerov.Fill(0.0);
	this->m_DisplacementField = DisplacementFieldType::New();
	this->m_DisplacementField->SetRegions(this->m_ReferenceSize);
	this->m_DisplacementField->SetSpacing(this->m_ReferenceSpacing);
	this->m_DisplacementField->SetOrigin(this->m_ReferenceOrigin);
	this->m_DisplacementField->SetDirection(this->m_ReferenceDirection);
	this->m_DisplacementField->Allocate();
	this->m_DisplacementField->FillBuffer(zerov);

	DisplacementType* dfbuffer = this->m_DisplacementField->GetBufferPointer();
	size_t nPix = this->m_DisplacementField->GetLargestPossibleRegion().GetNumberOfPixels();

	DisplacementType v;
	for( size_t c = 0; c < this->m_NumberOfTransforms; c++) {
		TransformComponentPointer tf = TransformComponentType::New();
		tf->SetOutputReference(this->m_DisplacementField);
		tf->SetCoefficientsVectorImage(this->m_Coefficients[c]);
		tf->Interpolate();
		TransformFieldConstPointer f = tf->GetField();

		const DisplacementType* cbuffer = f->GetBufferPointer();
		DisplacementType vc;
		for( size_t i = 0; i < nPix; i++) {
			vc = *(cbuffer + i);
			v = *(dfbuffer + i);
			*(dfbuffer + i)+= vc;
		}
	}
}

template< class TScalar, unsigned int NDimensions >
void
CompositeMatrixTransform<TScalar,NDimensions>
::SetPhysicalDomainInformation( const DomainBaseType* image ) {
	this->m_ReferenceSize = image->GetLargestPossibleRegion().GetSize();
	this->m_ReferenceSpacing = image->GetSpacing();
	this->m_ReferenceOrigin = image->GetOrigin();
	this->m_ReferenceDirection = image->GetDirection();
}

} // namespace rstk
