// --------------------------------------------------------------------------------------
// File:          rstkCoefficientsWriter.hxx
// Date:          Nov 27, 2014
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
