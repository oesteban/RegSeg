/* --------------------------------------------------------------------------------------
 * File:    GaussianModelDataEnergy.h
 * Date:    18/02/2012
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2012, Oscar Esteban - oesteban@die.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM) and
 Signal Processing Laboratory 5, EPFL (LTS5-EPFL).
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the names of the BIT-UPM and the LTS5-EPFL, nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef GAUSSIANMODELDATAENERGY_H_
#define GAUSSIANMODELDATAENERGY_H_

#include <vector>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkVariableSizeMatrix.h>


#include "DataEnergy.h"
#include "SparseDeformationField.h"
#include "SparseMultivariateInterpolator.h"
#include "IDWMultivariateInterpolator.h"
#include <itkNormalQuadEdgeMeshFilter.h>

namespace rstk
{

template< class TFixedImage, class TMovingSurface, class TEnergyValue = float >
class GaussianModelDataEnergy: public DataEnergy< TFixedImage, TMovingSurface, TEnergyValue > {
public:
	typedef GaussianModelDataEnergy                               Self;
	typedef DataEnergy<TFixedImage,TMovingSurface,TEnergyValue>   Superclass;
	typedef itk::SmartPointer<Self>                               Pointer;
	typedef itk::SmartPointer< const Self >                       ConstPointer;

	itkTypeMacro( GaussianModelDataEnergy, DataEnergy );
	itkNewMacro( Self );

	typedef TEnergyValue                                          EnergyValueType;

	typedef TFixedImage                                           FixedImageType;
	typedef typename FixedImageType::Pointer                      FixedImagePointer;
	typedef typename FixedImageType::ConstPointer                 FixedImageConstPointer;
	typedef typename FixedImageType::PixelType                    FeatureVectorType;
	itkStaticConstMacro(FeatureVectorSize, unsigned int, TFixedImage::PixelType::Dimension);
	itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::Dimension);

	typedef TMovingSurface                                        MovingSurfaceType;
	typedef typename MovingSurfaceType::Pointer                   MovingSurfacePointer;
	typedef typename MovingSurfaceType::ConstPointer              MovingSurfaceConstPointer;
	typedef typename MovingSurfaceType::PointType                 PointType;
	typedef std::vector< const MovingSurfaceType *>               MovingSurfaceContainer;

	typedef SparseDeformationField< MovingSurfaceType >           DisplacementFieldType;
	typedef typename DisplacementFieldType::Pointer               DisplacementFieldPointer;
	typedef typename DisplacementFieldType::ConstPointer          DisplacementFieldConstPointer;
	typedef std::vector< const DisplacementFieldType* >           DisplacementFieldContainer;
	typedef typename DisplacementFieldType::SparseVectorFieldType SparseVectorFieldType;
	typedef typename SparseVectorFieldType::Pointer               SparseVectorFieldPointer;
	typedef typename SparseVectorFieldType::ConstPointer          SparseVectorFieldConstPointer;
	typedef std::vector< const SparseVectorFieldType* >           SparseVectorFieldContainer;
	typedef typename SparseVectorFieldType::ParameterType         DisplacementVectorType;

	typedef SparseMultivariateInterpolator
			<MovingSurfaceType, FixedImageType>                   InterpolatorType;
	typedef typename InterpolatorType::Pointer                    InterpolatorPointer;
	typedef typename InterpolatorType::ConstPointer               InterpolatorConstPointer;

	typedef IDWMultivariateInterpolator
			<MovingSurfaceType, FixedImageType>                   DefaultInterpolatorType;

	typedef itk::Mesh< DisplacementVectorType, ImageDimension>    NormalsMeshType;
	typedef typename NormalsMeshType::Pointer                     NormalsMeshPointer;
	typedef typename NormalsMeshType::ConstPointer                NormalsMeshConstPointer;

	typedef itk::NormalQuadEdgeMeshFilter
			< MovingSurfaceType, NormalsMeshType >                NormalsCalculatorFilter;
	typedef typename NormalsCalculatorFilter::Pointer             NormalsCalculatorPointer;
	typedef typename NormalsCalculatorFilter::ConstPointer        NormalsCalculatorConstPointer;

	typedef itk::Array< double >                                  ParametersType;
    typedef itk::VariableSizeMatrix< double >                     MatrixType;

	//typedef itk::CovarianceSampleFilter< FixedImageType >

    enum SURF_TYPE {
    	INTERNAL,
    	EXTERNAL
    };

	itkSetMacro( ModelMean, FeatureVectorType );
	itkGetConstMacro( ModelMean, FeatureVectorType );
	itkSetMacro( ModelInvCovariance, MatrixType );
	itkGetConstMacro( ModelInvCovariance, MatrixType );


	itkSetConstObjectMacro( FixedImage, FixedImageType);
	itkGetConstObjectMacro( FixedImage, FixedImageType);

	itkSetConstObjectMacro( DisplacementField, DisplacementFieldType );
	itkGetConstObjectMacro( DisplacementField, DisplacementFieldType );

	EnergyValueType ComputeEnergy();
	EnergyValueType GetGradient();

	void SetModelParameters( FeatureVectorType mean, MatrixType cov );
protected:
	GaussianModelDataEnergy();
	virtual ~GaussianModelDataEnergy(){};

private:
	GaussianModelDataEnergy( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	SparseVectorFieldType* ComputeNormals();

	FixedImageConstPointer                    m_FixedImage;
	MovingSurfaceConstPointer                 m_MovingSurface;
	NormalsMeshConstPointer                   m_MovingNormals;
	DisplacementFieldContainer                m_DisplacementField;
	InterpolatorConstPointer                  m_Interpolator;

	FeatureVectorType                         m_ModelMean;
    MatrixType                                m_ModelInvCovariance;

}; // End of Class

} // End of namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "GaussianModelDataEnergy.txx"
#endif


#endif /* GAUSSIANMODELDATAENERGY_H_ */
