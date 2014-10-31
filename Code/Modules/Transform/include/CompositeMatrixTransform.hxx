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
void
CompositeMatrixTransform<TScalar,NDimensions>
::Compute() {
	if ((this->m_NumberOfTransforms == 0) || (this->m_Coefficients.size()!=this->m_NumberOfTransforms)) {
		itkExceptionMacro(<< "number of transforms is zero or it does not match the number of stored coefficients sets.");
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

	for( size_t c = 0; c < this->m_NumberOfTransforms; c++) {
		TransformComponentPointer tf = TransformComponentType::New();
		tf->SetPhysicalDomainInformation(this->m_DisplacementField);
		tf->SetCoefficientsImages(this->m_Coefficients[c]);
		tf->Interpolate();
		TransformFieldConstPointer f = tf->GetField();

		const typename TransformFieldType::PixelType* cbuffer = f->GetBufferPointer();

		for( size_t i = 0; i < nPix; i++) {
			*(dfbuffer + i) = *(dfbuffer + i) + *(cbuffer + i);
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
