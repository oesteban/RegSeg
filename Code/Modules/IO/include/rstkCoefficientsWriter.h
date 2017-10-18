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
