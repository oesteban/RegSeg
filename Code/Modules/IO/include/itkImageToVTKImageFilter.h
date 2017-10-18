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

#ifndef __itkImageToVTKImageFilter_h
#define __itkImageToVTKImageFilter_h

#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"

namespace itk
{

/** \class ImageToVTKImageFilter
 * \brief Converts an ITK image into a VTK image and plugs a
 *  itk data pipeline to a VTK datapipeline.
 *
 *  This class puts together an itkVTKImageExporter and a vtkImageImporter.
 *  It takes care of the details related to the connection of ITK and VTK
 *  pipelines. The User will perceive this filter as an adaptor to which
 *  an itk::Image can be plugged as input and a vtkImage is produced as
 *  output.
 *
 * \ingroup   ITKVtkGlue
 *
 * \wiki
 * \wikiexample{IO/ImageToVTKImageFilter,Display an ITK image}
 * \wikiexample{IO/itkVtkImageConvertDICOM,Uses a custom user matrix to align the image with DICOM physical space}
 * \endwiki
 */
template <typename TInputImage >
class ImageToVTKImageFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ImageToVTKImageFilter     Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToVTKImageFilter, ProcessObject);

  /** Some typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename InputImageType::ConstPointer  InputImagePointer;

  typedef VTKImageExport< InputImageType>        ExporterFilterType;
  typedef typename ExporterFilterType::Pointer   ExporterFilterPointer;

  /** Get the output in the form of a vtkImage.
      This call is delegated to the internal vtkImageImporter filter  */
  vtkImageData *  GetOutput() const;

  /** Set the input in the form of an itk::Image */
  using Superclass::SetInput;
  void SetInput( const InputImageType * );
  InputImageType * GetInput();

  /** Return the internal VTK image importer filter.
      This is intended to facilitate users the access
      to methods in the importer */
  vtkImageImport * GetImporter() const;

  /** Return the internal ITK image exporter filter.
      This is intended to facilitate users the access
      to methods in the exporter */
  ExporterFilterType * GetExporter() const;

  /** This call delegates the update to the importer */
  void Update();

protected:
  ImageToVTKImageFilter();
  virtual ~ImageToVTKImageFilter();

private:
  ImageToVTKImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);        //purposely not implemented

  ExporterFilterPointer       m_Exporter;
  vtkImageImport *            m_Importer;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToVTKImageFilter.hxx"
#endif

#endif
