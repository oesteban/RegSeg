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
#include <itkCommand.h>
#include <itkDataObjectDecorator.h>
#include <vector>       // std::vector
#include <iostream>     // std::cout


#include <jsoncpp/json/json.h>

#include "rstkMacro.h"
#include "ConfigurableObject.h"
#include "FunctionalBase.h"
#include "MeanFunctional.h"
#include "MahalanobisFunctional.h"
#include "OptimizerBase.h"
#include "SpectralGradientDescentOptimizer.h"
#include "SegmentationOptimizer.h"

#include "IterationJSONUpdate.h"
#include "IterationResultsWriterUpdate.h"


namespace bpo = boost::program_options;

namespace rstk {

template < typename TFixedImage, typename TTransform, typename TComputationalValue = double >
class ACWERegistrationMethod:
		public itk::ProcessObject,
		public ConfigurableObject
{
public:
	/** Standard class typedefs. */
	typedef ACWERegistrationMethod                            Self;
	typedef itk::ProcessObject                                Superclass;
	typedef itk::SmartPointer< Self >                         Pointer;
	typedef itk::SmartPointer< const Self >                   ConstPointer;
	typedef ConfigurableObject                                SettingsClass;
	typedef typename SettingsClass::SettingsMap               SettingsMap;
	typedef typename SettingsClass::SettingsDesc              SettingsDesc;
	typedef std::vector< SettingsMap >                        SettingsList;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );
	itkTypeMacro( ACWERegistrationMethod, itk::ProcessObject );

	typedef TComputationalValue                               InternalValuesType;

	typedef TFixedImage                                       ReferenceImageType;
	typedef typename ReferenceImageType::Pointer              ReferenceImagePointer;
	typedef typename ReferenceImageType::ConstPointer         ReferenceImageConstPointer;
	typedef typename ReferenceImageType::SizeType             GridSizeType;
	typedef std::vector< GridSizeType >                       GridSizeList;

	itkStaticConstMacro( Dimension, size_t, ReferenceImageType::ImageDimension );

	typedef itk::FixedArray< InternalValuesType, Dimension >  ValueArrayType;
	typedef itk::Vector< InternalValuesType, Dimension >      ValueVectorType;

	typedef TTransform                                        OutputTransformType;
	typedef typename OutputTransformType::Pointer             OutputTransformPointer;
	typedef itk::DataObjectDecorator<OutputTransformType>     DecoratedOutputTransformType;
	typedef typename DecoratedOutputTransformType::Pointer    DecoratedOutputTransformPointer;

	typedef FunctionalBase< ReferenceImageType >              FunctionalType;
	typedef typename FunctionalType::Pointer                  FunctionalPointer;
	typedef std::vector< FunctionalPointer >                  FunctionalList;
	typedef typename FunctionalType::ContourType              ContourType;
	typedef typename FunctionalType::ContourPointer           ContourPointer;
	typedef typename FunctionalType::ContourConstPointer      ContourConstPointer;
	typedef std::vector< ContourConstPointer >                PriorsList;


	typedef OptimizerBase< FunctionalType >                   OptimizerType;
	typedef typename OptimizerType::Pointer                   OptimizerPointer;
	typedef std::vector< OptimizerPointer >                   OptimizerList;
	typedef typename OptimizerType
			        ::InternalComputationValueType            OptCompValueType;
	typedef std::vector< OptCompValueType >                   OptCompValueList;

	typedef typename OptimizerType::SizeValueType             NumberValueType;
	typedef std::vector< NumberValueType >                    NumberValueList;

	//typedef MeanFunctional< ReferenceImageType >       DefaultFunctionalType;
	typedef SegmentationOptimizer< FunctionalType >    DefaultOptimizerType;
	typedef MahalanobisFunctional< ReferenceImageType >       DefaultFunctionalType;
	//typedef SpectralGradientDescentOptimizer
	//		                            < FunctionalType >    DefaultOptimizerType;

	typedef Json::Value                                       JSONRoot;
	typedef IterationJSONUpdate< OptimizerType >              JSONLoggerType;
	typedef typename JSONLoggerType::Pointer                  JSONLoggerPointer;


	typedef IterationResultWriterUpdate< OptimizerType > IterationWriterUpdate;
	typedef typename IterationWriterUpdate::Pointer      IterationWriterPointer;

	/** Codes of stopping conditions. */
	typedef enum {
		ALL_LEVELS_DONE,
		LEVEL_SETTING_ERROR,
		LEVEL_PROCESS_ERROR,
		INITIALIZATION_ERROR,
		OTHER_ERROR
	} StopConditionType;

	/** Stop condition return string type */
	typedef std::string                                             StopConditionReturnStringType;

	/** Stop condition internal string type */
	typedef std::ostringstream                                      StopConditionDescriptionType;

	itkSetInputMacro( FixedImage, ReferenceImageType );
	itkGetInputMacro( FixedImage, ReferenceImageType );

	//itkSetMacro( NumberOfLevels, size_t );
	void SetNumberOfLevels( size_t levels );
	itkGetConstMacro( NumberOfLevels, size_t );

	itkGetConstMacro( CurrentLevel, size_t );

	itkSetMacro( UseGridLevelsInitialization, bool );
	itkGetConstMacro( UseGridLevelsInitialization, bool );

	itkSetMacro( UseGridSizeInitialization, bool );
	itkGetConstMacro( UseGridSizeInitialization, bool );

	itkSetMacro( MaxGridSize, GridSizeType );
	itkGetConstMacro( MaxGridSize, GridSizeType );

	itkSetMacro( MinGridSize, GridSizeType );
	itkGetConstMacro( MinGridSize, GridSizeType );

	//itkGetConstMacro( StopConditionDescription, StopConditionDescriptionType );

	itkSetMacro( OutputPrefix, std::string );
	itkGetConstMacro( OutputPrefix, std::string );

	itkSetMacro( AutoSmoothing, bool );
	itkGetConstMacro( AutoSmoothing, bool );

	itkGetMacro( JSONRoot, JSONRoot);

	rstkSetVectorElement( GridSchedule, GridSizeType );
	rstkGetConstVectorElement( GridSchedule, GridSizeType );

	rstkVectorMethods( NumberOfIterations, NumberValueType );

	rstkVectorMethods( StepSize, OptCompValueType );
	rstkVectorMethods( Alpha, OptCompValueType );
	rstkVectorMethods( Beta, OptCompValueType );
	rstkVectorMethods( DescriptorRecomputationFreq, NumberValueType );

	rstkGetObjectList( Optimizer, OptimizerType );

	void AddShapePrior( const ContourType *prior ) { this->m_Priors.push_back( prior ); }

	// Methods inherited from the Configurable interface
	virtual void AddOptions( SettingsDesc& opts ) const {};

	void SetSettingsOfLevel( size_t l, SettingsMap& map );

	/** Returns the transform resulting from the registration process  */
	virtual const DecoratedOutputTransformType * GetOutput() const;


protected:
	ACWERegistrationMethod();
	~ACWERegistrationMethod() {}

	virtual void PrintSelf( std::ostream &os, itk::Indent indent ) const;
	virtual void GenerateData();

	void Initialize();
	void GenerateSchedule();
	void SetUpLevel( size_t level );
	void Stop( StopConditionType code, std::string msg );

	virtual void ParseSettings() {};
private:
	ACWERegistrationMethod( const Self & );
	void operator=( const Self & );

	size_t m_NumberOfLevels;
	size_t m_CurrentLevel;
	std::string m_OutputPrefix;
	bool m_UseGridLevelsInitialization;
	bool m_UseGridSizeInitialization;
	bool m_Initialized;
	bool m_AutoSmoothing;

	/* Common variables for optimization control and reporting */
	bool                          m_Stop;
	StopConditionType             m_StopCondition;
	StopConditionDescriptionType  m_StopConditionDescription;

	GridSizeList m_GridSchedule;
	GridSizeType m_MaxGridSize;
	GridSizeType m_MinGridSize;
	NumberValueList m_NumberOfIterations;

	FunctionalList m_Functionals;
	OptimizerList m_Optimizers;
	PriorsList m_Priors;
	SettingsList m_Config;
	OutputTransformPointer m_OutputTransform;
	OptCompValueList m_StepSize;
	OptCompValueList m_Alpha;
	OptCompValueList m_Beta;
	NumberValueList  m_DescriptorRecomputationFreq;

	JSONRoot m_JSONRoot;
	JSONLoggerPointer m_CurrentLogger;
	IterationWriterPointer m_ImageLogger;
};

} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ACWERegistrationMethod.hxx"
#endif

#endif /* ACWEREGISTRATIONMETHOD_H_ */
