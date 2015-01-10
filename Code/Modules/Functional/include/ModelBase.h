// --------------------------------------------------------------------------------------
// File:          ModelBase.h
// Date:          Dec 23, 2014
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
// This file is part of ACWEReg
//
// ACWEReg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACWEReg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ACWEReg.  If not, see <http://www.gnu.org/licenses/>.
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
	typedef typename PixelType::ValueType                                     PixelValueType;

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
	inline double Evaluate(const MeasurementVectorType & x, const RegionIdentifier roi) const {
		if (roi == this->m_NumberOfRegions - 1 ) return this->m_MaxEnergy;
		if (roi == this->m_NumberOfRegions - 2 ) return 0;

		return this->m_Memberships[roi]->Evaluate(x);
	}

    /** Set/Get priors
     *
     */
    virtual void SetPriorsMap(const PriorsImageType *priors) {
    	this->SetNthInput(1, const_cast<PriorsImageType *>(priors));
    }

    virtual const PriorsImageType * GetPriorsMap() {
    	return static_cast<const PriorsImageType*>(this->ProcessObject::GetInput(1));
    }

    /** Set/Get priors
     *
     */
    virtual void SetMask(const MaskType *mask) {
    	this->SetNthInput(2, const_cast<MaskType *>(mask));
    }

    virtual const MaskType * GetMask() {
    	return static_cast<const MaskType*>(this->ProcessObject::GetInput(2));
    }

	/** Get the length of the measurement vector */
	itkGetConstMacro(MeasurementVectorSize, MeasurementVectorSizeType);

	virtual MeasureTypeContainer GetRegionOffsetContainer() const = 0;
	virtual std::string PrintFormattedDescriptors() = 0;

	itkGetConstMacro(MaxEnergy, MeasureType);

protected:
	ModelBase();
	virtual ~ModelBase() {}

	virtual void PrintSelf(std::ostream & os, itk::Indent indent) const {
		Superclass::PrintSelf(os, indent);
		os << indent << "Length of measurement vectors: "
				<< m_MeasurementVectorSize << std::endl;
	}

	void GenerateInputRequestedRegion() {};
	virtual void InitializeMemberships() = 0;

	void SetMembership(const MembershipFunctionType* m, size_t i) {
		this->m_Memberships[i] = const_cast<MembershipFunctionType *>(m);
	}

	virtual MembershipFunctionType* GetNewFunction() = 0;

	inline MembershipFunctionType* GetFunction(size_t id) {
		return dynamic_cast<MembershipFunctionType*>( this->m_Memberships[id].GetPointer());
	}


	/** Method that construct the outputs */
	typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
	using Superclass::MakeOutput;
	DataObjectPointer MakeOutput(DataObjectPointerArraySizeType);

	MeasurementVectorSizeType m_MeasurementVectorSize;
	RegionIdentifier m_NumberOfRegions;
	MembershipFunctionsArray m_Memberships;
	MeasureType m_MaxEnergy;

private:
	ModelBase(const Self &);   //purposely not implemented
	void operator=(const Self &); //purposely not implemented
};

} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ModelBase.hxx"
#endif


#endif /* _MODELBASE_H_ */
