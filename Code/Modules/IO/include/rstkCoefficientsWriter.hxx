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

#ifndef RSTKCOEFFICIENTSWRITER_HXX_
#define RSTKCOEFFICIENTSWRITER_HXX_

#include "rstkCoefficientsWriter.h"
#include <fstream>

namespace rstk {

//
// Set the input mesh
//
template< typename TInputPointSet, typename TCoordRepType >
void
CoefficientsWriter<TInputPointSet, TCoordRepType>
::SetInput(const InputPointSetType *input)
{
  this->m_Input = input;
}

template< typename TInputPointSet, typename TCoordRepType >
void
CoefficientsWriter<TInputPointSet, TCoordRepType>
::SetCoefficientsImageArrayInput(const CoefficientsImageArray arr)
{
	typename InputPointSetType::Pointer input = InputPointSetType::New();

	CoeffImagePointer ref = arr[0];
	size_t npix = ref->GetLargestPossibleRegion().GetNumberOfPixels();

	PointType p;
	PixelType v;
	typename CoefficientsImageType::IndexType idx;

	for(size_t i = 0; i < npix; i++) {
		idx = ref->ComputeIndex(i);
		ref->TransformIndexToPhysicalPoint(idx, p);
		input->GetPoints()->InsertElement(i, p);

		for(size_t d = 0; d < Dimension; d++) {
			v[d] = arr[d]->GetPixel(idx);
		}
		input->GetPointData()->InsertElement(i, v);
	}

	this->m_Input = input;
}

//
// Write the input mesh to the output file
//
template< typename TInputPointSet, typename TCoordRepType >
void CoefficientsWriter<TInputPointSet, TCoordRepType>
::Update()
{
  this->GenerateData();
}

//
// Write the input mesh to the output file
//
template< typename TInputPointSet, typename TCoordRepType >
void CoefficientsWriter<TInputPointSet, TCoordRepType>
::Write()
{
  this->GenerateData();
}

template< typename TInputPointSet, typename TCoordRepType >
void CoefficientsWriter<TInputPointSet, TCoordRepType>::GenerateData() {
	if (this->m_FileName == "") {
		itkExceptionMacro("No FileName");
		return;
	}


	std::string::size_type idx;
	std::string extension = "";

	idx = this->m_FileName.rfind('.');

	if(idx != std::string::npos) {
	    extension = this->m_FileName.substr(idx+1);
	}

	if(extension.compare("vtu") == 0) {
		this->GenerateXMLData();
	} else {
		this->GenerateLegacyData();
	}
}

template< typename TInputPointSet, typename TCoordRepType >
void CoefficientsWriter<TInputPointSet, TCoordRepType>::GenerateLegacyData() {
	//
	// Write to output file
	//
	std::ofstream outputFile(this->m_FileName.c_str());


	if (!outputFile.is_open()) {
		itkExceptionMacro("Unable to open file\n"
				"outputFilename= " << this->m_FileName);
		return;
	}

	outputFile.imbue(std::locale::classic());
	outputFile << "# vtk DataFile Version 2.0" << std::endl;
	outputFile << "coefficients field" << std::endl;
	outputFile << "ASCII" << std::endl;
	outputFile << "FIELD parameters 2" << std::endl;

	// POINTS go first

	unsigned int numberOfPoints = this->m_Input->GetNumberOfPoints();
	outputFile << "LOCATIONS " << Dimension << " "
			<< numberOfPoints << " float" << std::endl;

	const PointsContainer *points = this->m_Input->GetPoints();

	if (points) {
		PointIterator pointIterator = points->Begin();
		PointIterator pointEnd = points->End();

		while (pointIterator != pointEnd) {
			PointType point = pointIterator.Value();

			outputFile << point[0] << " " << point[1];

			if (Dimension > 2) {
				outputFile << " " << point[2];
			} else {
				outputFile << " " << "0.0";
			}

			outputFile << std::endl;

			pointIterator++;
		}
	}
	outputFile << "DISPLACEMENTS " << Dimension << " "
			<< numberOfPoints << " float" << std::endl;

	const PointDataContainer *data = this->m_Input->GetPointData();
	if (data) {
		PointDataIterator pointIterator = data->Begin();
		PointDataIterator pointEnd = data->End();

		while (pointIterator != pointEnd) {
			PixelType point = pointIterator.Value();

			outputFile << point[0] << " " << point[1];

			if (Dimension > 2) {
				outputFile << " " << point[2];
			} else {
				outputFile << " " << "0.0";
			}

			outputFile << std::endl;
			pointIterator++;
		}
	}

	outputFile.close();
}

template< typename TInputPointSet, typename TCoordRepType >
void CoefficientsWriter<TInputPointSet, TCoordRepType>::GenerateXMLData() {
	//
	// Write to output file
	//
	std::ofstream outputFile(this->m_FileName.c_str());


	if (!outputFile.is_open()) {
		itkExceptionMacro("Unable to open file\n"
				"outputFilename= " << this->m_FileName);
		return;
	}

	size_t npoints = this->m_Input->GetNumberOfPoints();

	outputFile.imbue(std::locale::classic());
	outputFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">" << std::endl;
	outputFile << "\t<UnstructuredGrid>" << std::endl;
	outputFile << "\t\t<Piece NumberOfPoints=\""<< npoints << "\" NumberOfCells=\"0\">" << std::endl;

	outputFile << "\t\t\t<PointData Vectors=\"uk\">" << std::endl;
	outputFile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"uk\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

	const PointDataContainer *data = this->m_Input->GetPointData();
	if (data) {
		PointDataIterator pointIterator = data->Begin();
		PointDataIterator pointEnd = data->End();

		while (pointIterator != pointEnd) {
			PixelType point = pointIterator.Value();

			outputFile << "\t\t\t\t\t" << point[0] << " " << point[1];

			if (Dimension > 2) {
				outputFile << " " << point[2];
			} else {
				outputFile << " " << "0.0";
			}

			outputFile << std::endl;
			pointIterator++;
		}
	}
	outputFile << "\t\t\t\t</DataArray>" << std::endl;
	outputFile << "\t\t\t</PointData>" << std::endl;
	outputFile << "\t\t\t<CellData />" << std::endl;

	outputFile << "\t\t\t<Points>" << std::endl;
	outputFile << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

	// POINTS go first
	const PointsContainer *points = this->m_Input->GetPoints();
	if (points) {
		PointIterator pointIterator = points->Begin();
		PointIterator pointEnd = points->End();

		while (pointIterator != pointEnd) {
			PointType point = pointIterator.Value();

			outputFile << "\t\t\t\t\t" << point[0] << " " << point[1];

			if (Dimension > 2) {
				outputFile << " " << point[2];
			} else {
				outputFile << " " << "0.0";
			}

			outputFile << std::endl;
			pointIterator++;
		}
	}
	outputFile << "\t\t\t\t</DataArray>" << std::endl;
	outputFile << "\t\t\t</Points>" << std::endl;


	std::stringstream celltypes;
	std::stringstream cellconn;
	std::stringstream celloffsets;

	for(size_t i = 0; i < npoints; i++) {
		celltypes << "1 ";
		cellconn << i << " ";
		celloffsets << (i+1) << " ";
	}

	outputFile << "\t\t\t<Cells>" << std::endl;
	outputFile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
	outputFile << "\t\t\t\t\t" << celltypes.str() << std::endl;
	outputFile << "\t\t\t\t</DataArray>" << std::endl;

	outputFile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
	outputFile << "\t\t\t\t\t" << cellconn.str() << std::endl;
	outputFile << "\t\t\t\t</DataArray>" << std::endl;

	outputFile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
	outputFile << "\t\t\t\t\t" << celloffsets.str() << std::endl;
	outputFile << "\t\t\t\t</DataArray>" << std::endl;
	outputFile << "\t\t\t</Cells>" << std::endl;

	outputFile << "\t\t</Piece>" << std::endl;
	outputFile << "\t</UnstructuredGrid>" << std::endl;
	outputFile << "</VTKFile>" << std::endl;

	outputFile.close();
}

}  // namespace rstk

#endif /* RSTKCOEFFICIENTSWRITER_HXX_ */
