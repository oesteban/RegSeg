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

#ifndef ITERATIONRESULTSWRITERUPDATE_H_
#define ITERATIONRESULTSWRITERUPDATE_H_

#include <fstream>
#include <string>
#include <iostream>

#include "IterationUpdate.h"

#include <itkImageFileWriter.h>
#include <itkMeshFileWriter.h>

#include "DisplacementFieldFileWriter.h"
#include "rstkCoefficientsWriter.h"
#include "ComponentsFileWriter.h"
#include "DisplacementFieldComponentsFileWriter.h"


namespace rstk {

template< typename TOptimizer >
class IterationResultWriterUpdate: public IterationUpdate<TOptimizer>
{
public:
	typedef IterationResultWriterUpdate                 Self;
	typedef TOptimizer                                  OptimizerType;
	typedef typename itk::WeakPointer<OptimizerType>    OptimizerPointer;
	typedef IterationUpdate<TOptimizer>                 Superclass;
	typedef itk::SmartPointer<Self>                     Pointer;
	typedef itk::SmartPointer< const Self >             ConstPointer;

	typedef typename OptimizerType::FieldType                   FieldType;
	typedef typename OptimizerType::CoefficientsImageType       CoefficientsImageType;
	typedef typename OptimizerType::CoefficientsImageArray      CoefficientsImageArray;
	typedef typename OptimizerType::FunctionalType      		FunctionalType;
	typedef typename FunctionalType::ReferenceImageType         ReferenceImageType;
	typedef typename FunctionalType::ROIType                    ROIType;
	typedef typename FunctionalType::PriorsImageType            ProbabilityMapType;
	typedef typename FunctionalType::VectorContourType          VectorContourType;
	typedef typename OptimizerType::TransformType               TransformType;
	typedef typename TransformType::AltCoeffType                AltCoeffType;
	typedef itk::MeshFileWriter< VectorContourType >            ContourVectorWriterType;
	typedef typename ContourVectorWriterType::Pointer           ContourVectorWriterPointer;

	typedef rstk::DisplacementFieldFileWriter<FieldType> FieldWriter;
	typedef rstk::DisplacementFieldComponentsFileWriter<FieldType> ComponentsWriter;
	typedef typename itk::ImageFileWriter< CoefficientsImageType > CoefficientsWriter;
	typedef rstk::ComponentsFileWriter<ReferenceImageType>       ReferenceWriter;
	typedef rstk::ComponentsFileWriter<ProbabilityMapType>       MapWriter;

	itkTypeMacro( IterationResultWriterUpdate, IterationUpdate ); // Run-time type information (and related methods)
	itkNewMacro( Self );

	itkSetClampMacro( Verbosity, size_t, 0, 5 );
	itkGetConstMacro( Verbosity, size_t );

    itkSetMacro( Prefix, std::string );
    itkGetConstMacro( Prefix, std::string );

    void SetOptimizer( OptimizerType * optimizer ) {
      m_Optimizer = optimizer;
      m_Optimizer->AddObserver( itk::IterationEvent(), this );
      m_Optimizer->AddObserver( itk::StartEvent(), this );
      m_Optimizer->AddObserver( itk::EndEvent(), this );
      m_Optimizer->AddObserver( FunctionalModifiedEvent(), this );
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) {
    	if (this->m_Verbosity == 0 ) return;
    	std::stringstream ss;

    	if ( this->m_Prefix.size() > 0 ) {
    		char ch = *this->m_Prefix.rbegin();
    		if( ch != '_' )
    			this->m_Prefix.append( "_" );
    	}

		if( typeid( event ) == typeid( FunctionalModifiedEvent ) )  {
			if (this->m_Verbosity > 3 ) {
				ss.str("");
				ss << this->m_Prefix << "descriptors_" << std::setfill('0') << "lev" << this->m_Level << "_it"  << std::setw(3) << this->m_Optimizer->GetCurrentIteration() << ".json";
				std::string jsonstr = this->m_Optimizer->GetFunctional()->PrintFormattedDescriptors();

				std::ofstream outfile(ss.str().c_str());
				outfile << jsonstr;
				outfile.close();
			}
		}


    	if( typeid( event ) == typeid( itk::IterationEvent ) ) {

    		size_t nContours =this->m_Optimizer->GetFunctional()->GetCurrentContours().size();

    		if (this->m_Verbosity > 3 ) {
				ss.str("");
				typedef rstk::CoefficientsWriter< AltCoeffType > W;
				typename W::Pointer f = W::New();
				ss << this->m_Prefix << "uk_" << std::setfill('0') << "lev" << this->m_Level << "_it"  << std::setw(3) << this->m_Optimizer->GetCurrentIteration() << ".vtu";
				f->SetFileName( ss.str().c_str() );
				f->SetCoefficientsImageArrayInput(this->m_Optimizer->GetCoefficients());
				f->Update();

				ss.str("");
				typename W::Pointer fg = W::New();
				ss << this->m_Prefix << "gk_" << std::setfill('0') << "lev" << this->m_Level << "_it"  << std::setw(3) << this->m_Optimizer->GetCurrentIteration() << ".vtu";
				fg->SetFileName( ss.str().c_str() );
				fg->SetCoefficientsImageArrayInput(this->m_Optimizer->GetDerivativeCoefficients());
				fg->Update();
    		}

    		if( this->m_Verbosity > 1 ) {
				typename FunctionalType::VectorContourList grads = this->m_Optimizer->GetFunctional()->GetCurrentContours();
				for( size_t r = 0; r < nContours; r++ ) {
					ContourVectorWriterPointer wc = ContourVectorWriterType::New();
					std::stringstream ss;
					ss << this->m_Prefix << "gi_lev" << this->m_Level << "_it" << std::setfill('0')<< std::setw(3) << this->m_Optimizer->GetCurrentIteration() << std::setw(2) << "_cont"<< r << ".vtk";
					wc->SetFileName( ss.str().c_str() );
					wc->SetFileTypeAsASCII();
					wc->SetInput( grads[r] );
					wc->Update();
				}
       		}

    		if ( this->m_Verbosity > 2 ) {
				ss.str("");
				ss << this->m_Prefix << "regions_lev" << this->m_Level << "_it" << std::setfill('0')<<std::setw(3) << this->m_Optimizer->GetCurrentIteration() << ".nii.gz";
				typename MapWriter::Pointer wr = MapWriter::New();
				wr->SetInput( this->m_Optimizer->GetFunctional()->GetCurrentMaps());
				wr->SetFileName(ss.str().c_str() );
				wr->Update();

    		}

    	}

    	if (typeid( event ) == typeid( itk::StartEvent )) {
    		if ( this->m_Verbosity > 0 ) {
				typedef itk::ImageFileWriter< ROIType > WriteROI;
				typename WriteROI::Pointer w = WriteROI::New();
				std::stringstream ss;
				ss << this->m_Prefix << "initial_seg_" << this->m_Level << ".nii.gz";
				w->SetFileName( ss.str().c_str() );
				w->SetInput( this->m_Optimizer->GetFunctional()->GetCurrentRegions() );
				w->Update();
    		}

    		if ( this->m_Verbosity > 1 ) {
				typename ReferenceWriter::Pointer w2 = ReferenceWriter::New();
				std::stringstream ss;
				ss << this->m_Prefix << "reference_lev" << this->m_Level << ".nii.gz";
				w2->SetFileName( ss.str().c_str() );
				w2->SetInput( this->m_Optimizer->GetFunctional()->GetReferenceImage() );
				w2->Update();
    		}
    	}
    }

protected:
    IterationResultWriterUpdate(): m_Verbosity(1), m_Prefix("") {}
    ~IterationResultWriterUpdate(){}

private:
    IterationResultWriterUpdate( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemented

	OptimizerPointer   m_Optimizer;
	std::string        m_Prefix;
	size_t             m_Verbosity;
};

} // end namespace rstk

#endif /* ITERATIONRESULTSWRITERUPDATE_H_ */
