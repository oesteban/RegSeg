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

#ifndef _MODELBASE_H_
#define _MODELBASE_H_

#include <itkImageTransformer.h>
#include <itkMembershipFunctionBase.h>
#include <itkMeasurementVectorTraits.h>
#include <itkNumericTraitsCovariantVectorPixel.h>
#include <itkVectorImageToImageAdaptor.h>

namespace rstk {

template< typename TInputVectorImage, typename TPriorsPrecisionType = float>
class ModelBase: public itk::ImageTransformer< TInputVectorImage > {

public:
	/** Standard class typedefs */
	typedef ModelBase                                  Self;
	typedef itk::ImageTransformer< TInputVectorImage > Superclass;
	typedef itk::SmartPointer< Self >                  Pointer;
	typedef itk::SmartPointer< const Self >            ConstPointer;

	/** Standard macros */
	itkTypeMacro(ModelBase, itk::ImageTransformer);
	itkStaticConstMacro(Dimension, unsigned int, TInputVectorImage::ImageDimension);

	typedef TInputVectorImage                                                 InputImageType;
	typedef typename InputImageType::Pointer                                  InputImagePointer;
	typedef typename InputImageType::SizeType                                 SizeType;
	typedef typename InputImageType::SpacingType                              SpacingType;
	typedef typename InputImageType::IndexType                                IndexType;
	typedef typename InputImageType::PixelType                                PixelType;
	typedef typename InputImageType::PointType                                PointType;
	typedef typename InputImageType::RegionType                               RegionType;
	typedef typename InputImageType::InternalPixelType                        PixelValueType;

	/** MeasurementVector typedef support */
	typedef double MeasureType;
	typedef PixelType                                                         MeasurementVectorType;
	typedef typename itk::Array<MeasureType>                                  MeasureTypeContainer;

	/** Typedef for the length of each measurement vector */
	typedef unsigned int MeasurementVectorSizeType;
	typedef unsigned int RegionIdentifier;

	typedef TPriorsPrecisionType                                              PriorsPrecisionType;
	typedef itk::VectorImage< PriorsPrecisionType, Dimension >                PriorsImageType;
	typedef typename PriorsImageType::PixelType                               PriorsPixelType;
	typedef typename PriorsImageType::Pointer                                 PriorsImagePointer;
	typedef itk::ImageRegionConstIterator< PriorsImageType >                  PriorsImageIteratorType;

	typedef itk::VectorImageToImageAdaptor<PriorsPrecisionType, Dimension>    PriorsAdaptor;
	typedef typename PriorsAdaptor::Pointer                                   PriorsAdaptorPointer;

	typedef itk::Image< PriorsPrecisionType, Dimension >                      MaskType;
	typedef typename MaskType::Pointer                                        MaskPointer;

	typedef itk::Statistics::MembershipFunctionBase< MeasurementVectorType >  MembershipFunctionType;
	typedef typename MembershipFunctionType::Pointer                          MembershipFunctionPointer;
	typedef typename MembershipFunctionType::ConstPointer                     MembershipFunctionConstPointer;
	typedef std::vector< MembershipFunctionPointer >                          MembershipFunctionsArray;
	typedef itk::SimpleDataObjectDecorator< MembershipFunctionsArray >        MembershipFunctionsObject;

	typedef typename Superclass::DataObjectPointer                            DataObjectPointer;
	typedef std::vector< DataObjectPointer >                                  DataObjectPointerArray;

	/** Method to get membership score (discriminant score) of an entity
	 * or measurement. Evaluate() maps from a vector measurement type
	 * to a real number. */
	virtual double Evaluate(const MeasurementVectorType & x, const RegionIdentifier roi) const = 0;

    /** Set/Get priors
     *
     */
    virtual void SetPriorsMap(const PriorsImageType *priors) {
    	this->SetNthInput(1, const_cast<PriorsImageType *>(priors));
    }

    virtual const PriorsImageType * GetPriorsMap() {
    	return static_cast<const PriorsImageType*>(this->Superclass::ProcessObject::GetInput(1));
    }

    /** Set/Get priors
     *
     */
    virtual void SetMask(const MaskType *mask) {
    	this->SetNthInput(2, const_cast<MaskType *>(mask));
    }

    virtual const MaskType * GetMask() {
    	return static_cast<const MaskType*>(this->Superclass::ProcessObject::GetInput(2));
    }

	/** Get the length of the measurement vector */
	itkGetConstMacro(MeasurementVectorSize, MeasurementVectorSizeType);

	virtual MeasureTypeContainer GetRegionOffsetContainer() const = 0;
	virtual std::string PrintFormattedDescriptors() = 0;
	virtual void ReadDescriptorsFromFile(std::string filename) = 0;

	itkGetConstMacro(MaxEnergy, MeasureType);

	itkSetMacro(NumberOfSpecialRegions, size_t);
	itkGetMacro(NumberOfSpecialRegions, size_t);

protected:
	ModelBase();
	virtual ~ModelBase() {}

	virtual void PrintSelf(std::ostream & os, itk::Indent indent) const override {
		Superclass::PrintSelf(os, indent);
		os << indent << "Length of measurement vectors: "
				<< m_MeasurementVectorSize << std::endl;
	}

	void GenerateInputRequestedRegion() override {};
	virtual void InitializeMemberships() = 0;

	void SetMembership(const MembershipFunctionType* m, size_t i) {
		this->m_Memberships[i] = const_cast<MembershipFunctionType *>(m);
	}

	virtual MembershipFunctionType* GetNewFunction() = 0;

	inline MembershipFunctionType* GetFunction(size_t id) {
		return dynamic_cast<MembershipFunctionType*>( this->m_Memberships[id].GetPointer());
	}


	/** Method that construct the outputs */
	typedef itk::ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
	using Superclass::MakeOutput;
	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) override;

	MeasurementVectorSizeType   m_MeasurementVectorSize;
	RegionIdentifier            m_NumberOfRegions;
	MembershipFunctionsArray    m_Memberships;
	MeasureType                 m_MaxEnergy;
	size_t                      m_NumberOfSpecialRegions;

private:
	ModelBase(const Self &);   //purposely not implemented
	void operator=(const Self &); //purposely not implemented
};

} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ModelBase.hxx"
#endif


#endif /* _MODELBASE_H_ */
