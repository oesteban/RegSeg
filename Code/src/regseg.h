/* --------------------------------------------------------------------------------------
 * File:    regseg.h
 * Date:    02/10/2012
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

#ifndef REGSEG_H_
#define REGSEG_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMesh.h>
#include <itkMeshFileReader.h>
#include <itkImageToVectorImageFilter.h>


namespace bpo = boost::program_options;

typedef float       PixelType;
const unsigned int Dimension = 3;
typedef itk::VectorImage< PixelType, Dimension >    ImageVectorType;
typedef itk::ImageFileReader< ImageVectorType >     ImageReader;
typedef itk::ImageFileWriter< ImageVectorType >     ImageWriter;

typedef itk::Image< PixelType, Dimension >          ImageComponentType;
typedef itk::ImageFileReader<ImageComponentType>    ImageComponentReader;
typedef itk::ImageToVectorImageFilter
		          < ImageComponentType >            ImageToVectorFilterType;

typedef itk::Mesh< PixelType, Dimension >           SurfaceType;
typedef itk::MeshFileReader< SurfaceType >          SurfaceReader;
typedef std::vector< const SurfaceType* >           SurfaceVector;


int main(int argc, char *argv[]);

#endif /* REGSEG_H_ */
