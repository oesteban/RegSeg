/* --------------------------------------------------------------------------------------
 * File:    DataEnergy.h
 * Date:    19/02/2012
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 */

// Copyright (c) 2012, Oscar Esteban - oesteban@die.upm.es
// with Biomedical Image Technology, UPM (BIT-UPM) and
// Signal Processing Laboratory 5, EPFL (LTS5-EPFL).
// All rights reserved.
//
// This file is part of ACWE-Reg
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef DATAENERGY_H_
#define DATAENERGY_H_

#include "EnergyTerm.h"
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkMesh.h>

namespace rstk
{

template< class TFixedImage, class TMovingSurface, class TDeformationField, class TEnergyValue = float >
class DataEnergy:public EnergyTerm< TEnergyValue > {
public:
	typedef DataEnergy                                         Self;
	typedef EnergyTerm < TEnergyValue >                        Superclass;
	typedef itk::SmartPointer<Self>                            Pointer;
	typedef itk::SmartPointer< const Self >                    ConstPointer;

	itkTypeMacro( DataEnergy, EnergyTerm );

	typedef TEnergyValue                                       EnergyValueType;

	typedef TFixedImage                                        FixedImageType;
	typedef typename FixedImageType::Pointer                   FixedImagePointer;
	typedef typename FixedImageType::ConstPointer              FixedImageConstPointer;
	typedef typename FixedImageType::PixelType                 MeasurementVectorType;

	typedef TMovingSurface                                     MovingSurfaceType;
	typedef typename MovingSurfaceType::Pointer                MovingSurfacePointer;
	typedef typename MovingSurfaceType::ConstPointer           MovingSurfaceConstPointer;
	typedef std::vector< const MovingSurfaceType *>            MovingSurfaceContainer;


	itkSetConstObjectMacro( FixedImage, FixedImageType);
	itkGetObjectMacro( FixedImage, FixedImageType);

	virtual EnergyValueType ComputeEnergy() = 0;
	virtual float GetGradient() = 0;

protected:
	DataEnergy();
	virtual ~DataEnergy(){};

	FixedImageConstPointer m_FixedImage;
private:


	DataEnergy( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

}; // End of Class

} // End of namespace rstk


#endif /* DATAENERGY_H_ */
