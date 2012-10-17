// --------------------------------------------------------------------------
// File:             EnergyOptimizerBase.h
// Date:             16/10/2012
// Author:           code@oscaresteban.es (Oscar Esteban, OE)
// Version:          0.1
// License:          BSD
// --------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
// 
// This file is part of ACWE-Registration
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the names of the LTS5-EFPL and the BIT-UPM, nor the names of its
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

#ifndef ENERGYOPTIMIZERBASE_H_
#define ENERGYOPTIMIZERBASE_H_

#include <itkObject.h>
#include <itkOptimizerParameters.h>
#include "EnergyBase.h"
#include <itkIntTypes.h>

namespace rstk {
/** \class EnergyOptimizerBase
 * \brief Abstract base for energy minimizers.
 */
class EnergyOptimizerBase : public itk::Object {
public:
  /** Standard class typedefs. */
  typedef EnergyOptimizerBase                              Self;
  typedef Object                                           Superclass;
  typedef itk::SmartPointer< Self >                        Pointer;
  typedef itk::SmartPointer< const Self >                  ConstPointer;

  itkTypeMacro(EnergyOptimizerBase, itk::Object );

  typedef itk::OptimizerParameters< double >               ScalesType;     //  Scale type.
  typedef itk::OptimizerParameters< double >               ParametersType; //  Parameters type.

  /** Metric function type */
  typedef EnergyBase                                       EnergyType;
  typedef EnergyType::Pointer                              EnergyPointer;

  /** Number of parameters type */
  typedef EnergyType::NumberOfParametersType               NumberOfParametersType;

  /** Measure type */
  typedef EnergyType::MeasureType                          MeasureType;

  /** Internal computation value type */
  typedef EnergyType::InternalComputationValueType         InternalComputationValueType;

  /** Accessors for Energy */
  itkGetObjectMacro( Energy, EnergyType );
  itkSetObjectMacro( Energy, EnergyType );

  /** Accessor for Energy value. Returns the value
   *  stored in m_CurrentEnergyValue from the most recent
   *  call to evaluate the Energy. */
  itkGetConstReferenceMacro( CurrentEnergyValue, MeasureType );

  /** Set current parameters scaling. */
  itkSetMacro( Scales, ScalesType );

  /** Get current parameters scaling. */
  itkGetConstReferenceMacro( Scales, ScalesType );

  /** Get whether scales are identity. Cannot be set */
  itkGetConstReferenceMacro( ScalesAreIdentity, bool );

  /** Set the number of threads to use when threading.
   * The default is the global default number of threads
   * returned from itkMultiThreader. */
  virtual void SetNumberOfThreads( itk::ThreadIdType number );

  /** Get the number of threads set to be used. */
  itkGetConstReferenceMacro( NumberOfThreads, itk::ThreadIdType );

  /** Get a reference to the current position of the optimization.
   * This returns the parameters from the assigned Energy, since the optimizer
   * itself does not store a position. */
  const ParametersType & GetCurrentPosition();

  /** Run the optimization.
   * \note Derived classes must override and call this superclass method, then
   * perform any additional initialization before performing optimization. */
  virtual void Start();

protected:

  /** Default constructor */
  EnergyOptimizerBase();
  virtual ~EnergyOptimizerBase();

  EnergyPointer                 m_Energy;
  itk::ThreadIdType             m_NumberOfThreads;

  /** Energy measure value at a given iteration, as most recently evaluated. */
  MeasureType                   m_CurrentEnergyValue;

  /** Scales. Size is expected to be == Energy->GetNumberOfLocalParameters().
   * See the main documentation for more details. */
  ScalesType                    m_Scales;

  /** Flag to avoid unnecessary arithmetic when scales are identity. */
  bool                          m_ScalesAreIdentity;

  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  EnergyOptimizerBase( const Self & );  //purposely not implemented
  void operator=( const Self& );  //purposely not implemented
};

} // end namespace rstk


#endif /* ENERGYOPTIMIZERBASE_H_ */
