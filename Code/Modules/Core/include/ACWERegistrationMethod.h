// --------------------------------------------------------------------------------------
// File:          ACWERegistrationMethod.h
// Date:          Jan 29, 2014
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

#ifndef ACWEREGISTRATIONMETHOD_H_
#define ACWEREGISTRATIONMETHOD_H_

#include <itkProcessObject.h>

#include "FunctionalBase.h"
#include "SpectralOptimizer.h"

namespace rstk {

template < typename TFixedImage, typename TTransform >
class ACWERegistrationMethod: public itk::ProcessObject
{
public:
	/** Standard class typedefs. */
	typedef ACWERegistrationMethod                            Self;
	typedef itk::ProcessObject                                Superclass;
	typedef itk::SmartPointer< Self >                         Pointer;
	typedef itk::SmartPointer< const Self >                   ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );
	itkTypeMacro( ACWERegistrationMethod, itk::ProcessObject );

	typedef TFixedImage                                       ReferenceImageType;
	typedef typename ReferenceImageType::Pointer              ReferenceImagePointer;

	typedef FunctionalBase< ReferenceImageType >              FunctionalType;
	typedef typename FunctionalType::Pointer                  FunctionalPointer;

	typedef SpectralOptimizer< FunctionalType >               OptimizerType;
	typedef typename OptimizerType::Pointer                   OptimizerPointer;

	itkSetMacro( NumberOfLevels, size_t );
	itkGetConstMacro( NumberOfLevels, size_t );

protected:
	ACWERegistrationMethod();
	~ACWERegistrationMethod() {}

	virtual void GenerateData();

	void Initialize();

private:
	ACWERegistrationMethod( const Self & );
	void operator=( const Self & );

	size_t m_NumberOfLevels;
};

} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ACWERegistrationMethod.hxx"
#endif

#endif /* ACWEREGISTRATIONMETHOD_H_ */
