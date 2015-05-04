// --------------------------------------------------------------------------------------
// File:          rstkCoefficientsWriter.h
// Date:          Nov 27, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2014, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of RegSeg
//
// RegSeg is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RegSeg is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RegSeg.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef RSTKCOEFFICIENTSWRITER_H_
#define RSTKCOEFFICIENTSWRITER_H_

#include <itkObject.h>
#include <itkPointSet.h>

namespace rstk {
/** \class rstkCoefficientsWriter
 * \brief
 * Modifies the standard VTKPolyDataWriter behavior to write files
 * compatible with freeview (from FreeSurfer).
 * \ingroup ITKMesh
 */

template<class TInputPointSet, typename TCoordRepType = float>
class CoefficientsWriter: public itk::Object {
public:
	/** Standard "Self" typedef. */
	typedef CoefficientsWriter Self;
	typedef itk::Object Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Method for creation through the object factory */
	itkNewMacro (Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(CoefficientsWriter, itk::Object);

	itkStaticConstMacro( Dimension, unsigned int, TInputPointSet::PointDimension );

	typedef TCoordRepType ScalarType;

	/** Hold on to the type information specified by the template parameters.
	 */
	typedef TInputPointSet InputPointSetType;
	typedef typename InputPointSetType::PixelType PixelType;
	typedef typename InputPointSetType::PointType PointType;
	typedef typename InputPointSetType::PointIdentifier PointIdentifier;

	/** Some convenient typedefs. */
	typedef typename InputPointSetType::ConstPointer InputMeshPointer;

	typedef typename InputPointSetType::PointsContainer PointsContainer;
	typedef typename InputPointSetType::PointDataContainer PointDataContainer;

	typedef typename PointsContainer::ConstIterator PointIterator;
	typedef typename PointDataContainer::ConstIterator PointDataIterator;

    typedef itk::Image< ScalarType, Dimension >                                 CoefficientsImageType;
    typedef typename CoefficientsImageType::Pointer                             CoeffImagePointer;
    typedef typename CoefficientsImageType::ConstPointer                        CoeffImageConstPointer;
    typedef itk::FixedArray< CoeffImagePointer, Dimension >                     CoefficientsImageArray;

	void Update();
	void Write();
	void SetInput(const InputPointSetType *input);

	void SetCoefficientsImageArrayInput(const CoefficientsImageArray arr);

	/** Set/Get the name of the file where data are written. */
	itkSetStringMacro(FileName);
	itkGetStringMacro(FileName);

protected:
	CoefficientsWriter() :
		m_Input(NULL),
		m_FileName(""),
		Superclass() {}
	virtual ~CoefficientsWriter() {}

	virtual void GenerateData();
	virtual void GenerateLegacyData();
	virtual void GenerateXMLData();

	std::string m_FileName;
	InputMeshPointer m_Input;

	void PrintSelf(std::ostream & os, itk::Indent indent) const {
		Superclass::PrintSelf(os, indent);
		os << indent << "FileName: " << this->m_FileName << std::endl;
	}


private:
	CoefficientsWriter(const Self &); //purposely not implemented
	void operator=(const Self &);    //purposely not implemented
};
}  // namespace rstk

#ifndef ITK_MANUAL_INSTANTIATION
#include "rstkCoefficientsWriter.hxx"
#endif

#endif /* RSTKCOEFFICIENTSWRITER_H_ */
