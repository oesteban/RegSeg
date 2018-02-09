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
#include "OptimizerBase.h"
#include "SpectralGradientDescentOptimizer.h"
#include "SegmentationOptimizer.h"
#include "CompositeMatrixTransform.h"

#include "IterationJSONUpdate.h"
#include "IterationStdOutUpdate.h"
#include "IterationResultsWriterUpdate.h"

namespace bpo = boost::program_options;

namespace rstk {

template < typename TFixedImage, typename TTransform, typename TComputationalValue = float >
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

	typedef TTransform                                        TransformType;
	typedef typename TransformType::Pointer                   TransformPointer;
	typedef std::vector< TransformPointer >                   TransformList;

	typedef rstk::CompositeMatrixTransform
			            < InternalValuesType, Dimension >     OutputTransformType;
	typedef typename OutputTransformType::Pointer             OutputTransformPointer;
	typedef itk::DataObjectDecorator<OutputTransformType>     DecoratedOutputTransformType;
	typedef typename DecoratedOutputTransformType::Pointer    DecoratedOutputTransformPointer;
	typedef typename OutputTransformType::
			                            DisplacementFieldType OutputFieldType;
	typedef typename OutputTransformType::
			                         DisplacementFieldPointer OutputFieldPointer;
	typedef typename OutputTransformType::OutputVectorType    OutputVectorType;

	typedef FunctionalBase< ReferenceImageType >              FunctionalType;
	typedef typename FunctionalType::Pointer                  FunctionalPointer;
	typedef std::vector< FunctionalPointer >                  FunctionalList;
	typedef typename FunctionalType::ScalarContourCopyType    ContourCopyType;
	typedef typename ContourCopyType::Pointer                 ContourCopyPointer;
	typedef typename FunctionalType::VectorContourCopyType    ShapeCopyType;
	typedef typename ShapeCopyType::Pointer                   ShapeCopyPointer;
	typedef typename FunctionalType::Vector2ScalarCopyType    Shape2PriorCopyType;
	typedef typename Shape2PriorCopyType::Pointer             Shape2PriorCopyPointer;

	typedef typename FunctionalType::ROIType                  ROIType;
	typedef typename FunctionalType::ScalarContourType        PriorsType;
	typedef typename PriorsType::Pointer                      PriorPointer;
	typedef typename PriorsType::ConstPointer                 PriorConstPointer;
	typedef std::vector< PriorConstPointer >                  PriorsList;

	typedef typename FunctionalType::VectorContourType        VectorContourType;
	typedef typename VectorContourType::Pointer               ShapePointer;
	typedef typename VectorContourType::ConstPointer          ShapeConstPointer;
	typedef std::vector< ShapeConstPointer >                  ShapesList;

	typedef typename FunctionalType::ProbabilityMapType       FixedMaskType;
	typedef typename FixedMaskType::ConstPointer              FixedMaskConstPointer;


	typedef OptimizerBase< FunctionalType >                   OptimizerType;
	typedef typename OptimizerType::Pointer                   OptimizerPointer;
	typedef std::vector< OptimizerPointer >                   OptimizerList;
	typedef typename OptimizerType
			        ::InternalComputationValueType            OptCompValueType;
	typedef std::vector< OptCompValueType >                   OptCompValueList;

	typedef typename OptimizerType::SizeValueType             NumberValueType;
	typedef std::vector< NumberValueType >                    NumberValueList;

	typedef typename OptimizerType::FieldType                 FieldType;
	typedef typename FieldType::Pointer                       FieldPointer;
	typedef typename FieldType::ConstPointer                  FieldConstPointer;
	typedef std::vector< FieldConstPointer >                  FieldList;

	typedef SpectralGradientDescentOptimizer
			                            < FunctionalType >    DefaultOptimizerType;

	typedef Json::Value                                       JSONRoot;
	typedef IterationJSONUpdate< OptimizerType >              JSONLoggerType;
	typedef typename JSONLoggerType::Pointer                  JSONLoggerPointer;

	typedef IterationStdOutUpdate< OptimizerType >            STDOutLoggerType;
	typedef typename STDOutLoggerType::Pointer                STDOutLoggerPointer;

	typedef IterationResultWriterUpdate< OptimizerType >      IterationWriterUpdate;
	typedef typename IterationWriterUpdate::Pointer           IterationWriterPointer;

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

	itkSetConstObjectMacro( FixedMask, FixedMaskType );
	itkGetConstObjectMacro( FixedMask, FixedMaskType );

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

	void SetReferenceNames(const std::vector< std::string > s) { this->m_ReferenceNames = std::vector<std::string>(s); }
	void SetPriorsNames(const std::vector< std::string > s) { this->m_PriorsNames = std::vector<std::string>(s); }
	void SetTargetNames(const std::vector< std::string > s) { this->m_TargetNames = std::vector<std::string>(s); }

	itkSetMacro( OutputPrefix, std::string );
	itkGetConstMacro( OutputPrefix, std::string );

	itkSetMacro( AutoSmoothing, bool );
	itkGetConstMacro( AutoSmoothing, bool );

	itkSetClampMacro( Verbosity, size_t, 0, 5 );
	itkGetConstMacro( Verbosity, size_t );

	itkSetClampMacro( TransformNumberOfThreads, size_t, 1, ITK_MAX_THREADS );
	itkGetConstMacro( TransformNumberOfThreads, size_t );

	itkGetConstObjectMacro(OutputTransform, OutputTransformType );
	itkGetConstObjectMacro(OutputInverseTransform, OutputTransformType );
	itkGetConstObjectMacro(DisplacementField, FieldType );
	itkGetConstObjectMacro(InverseDisplacementField, FieldType );

	// itkGetConstMacro(Transforms, TransformList);
	itkGetMacro( JSONRoot, JSONRoot);

	rstkSetVectorElement( GridSchedule, GridSizeType );
	rstkGetConstVectorElement( GridSchedule, GridSizeType );

	rstkVectorMethods( NumberOfIterations, NumberValueType );

	rstkVectorMethods( StepSize, OptCompValueType );
	rstkVectorMethods( Alpha, OptCompValueType );
	rstkVectorMethods( Beta, OptCompValueType );
	rstkVectorMethods( DescriptorRecomputationFreq, NumberValueType );

	itkGetObjectMacro(Optimizer, OptimizerType);

	// rstkGetObjectListWithLast( Transform, TransformType );
	// rstkGetObjectListWithLast( Optimizer, OptimizerType );
	// rstkGetObjectListWithLast( Functional, FunctionalType );

	void AddShapeTarget( const PriorsType *surf ) { this->m_Target.push_back( surf ); }

	// Methods inherited from the Configurable interface
	virtual void AddOptions( SettingsDesc& opts ) const {};

	void SetSettingsOfLevel( size_t l, SettingsMap& map );

	/** Returns the transform resulting from the registration process  */
	virtual const DecoratedOutputTransformType * GetOutput() const;

	const FieldType* GetCurrentDisplacementField() const {
		return static_cast<const FieldType* >(this->m_Optimizer->GetCurrentDisplacementField());
	}

	FieldList GetCoefficientsField();

	PriorsList GetCurrentContours() const { return m_CurrentContours; }

	const ROIType* GetCurrentRegion( size_t contour_id ) const {
		return this->m_Functional->GetCurrentRegion( contour_id );
	}

protected:
	ACWERegistrationMethod();
	~ACWERegistrationMethod() {}

	virtual void PrintSelf( std::ostream &os, itk::Indent indent ) const override;
	virtual void GenerateData() override;

	void Initialize();
	void GenerateSchedule();
	void GenerateFinalDisplacementField();
	void ConcatenateFields( size_t level = 0 );
	void SetUpLevel( size_t level );
	void Stop( StopConditionType code, std::string msg );

	virtual void ParseSettings() override {};
private:
	ACWERegistrationMethod( const Self & );
	void operator=( const Self & );

	FixedMaskConstPointer m_FixedMask;

	size_t m_NumberOfLevels;
	size_t m_CurrentLevel;
	std::string m_OutputPrefix;
	bool m_UseGridLevelsInitialization;
	bool m_UseGridSizeInitialization;
	bool m_UseCustomGridSize;
	bool m_Initialized;
	bool m_AutoSmoothing;

	/* Common variables for optimization control and reporting */
	bool                          m_Stop;
	StopConditionType             m_StopCondition;
	StopConditionDescriptionType  m_StopConditionDescription;

	GridSizeList m_GridSchedule;
	GridSizeList m_FactorsSchedule;
	GridSizeType m_MaxGridSize;
	GridSizeType m_MinGridSize;
	NumberValueList m_NumberOfIterations;

	// TransformList m_Transforms;
	FunctionalPointer m_Functional;
	OptimizerPointer m_Optimizer;

	// FunctionalList m_Functionals;
	// OptimizerList m_Optimizers;
	PriorsList m_Target;
	PriorsList m_CurrentContours;
	SettingsList m_Config;
	OutputTransformPointer m_OutputTransform;
	OutputTransformPointer m_OutputInverseTransform;
	FieldPointer m_DisplacementField;
	FieldPointer m_InverseDisplacementField;
	FieldList m_CoefficientsContainer;
	OptCompValueList m_StepSize;
	OptCompValueList m_Alpha;
	OptCompValueList m_Beta;
	NumberValueList  m_DescriptorRecomputationFreq;

	JSONRoot m_JSONRoot;
	JSONLoggerPointer m_CurrentLogger;
	IterationWriterPointer m_ImageLogger;
	STDOutLoggerPointer m_OutLogger;

	size_t m_Verbosity;

	size_t m_TransformNumberOfThreads;

	std::vector< std::string > m_ReferenceNames;
	std::vector< std::string > m_PriorsNames;
	std::vector< std::string > m_TargetNames;
};

} // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ACWERegistrationMethod.hxx"
#endif

#endif /* ACWEREGISTRATIONMETHOD_H_ */
