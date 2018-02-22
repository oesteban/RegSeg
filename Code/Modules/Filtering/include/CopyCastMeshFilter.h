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

#ifndef SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_COPYCASTMESHFILTER_H_
#define SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_COPYCASTMESHFILTER_H_

#include <itkMeshToMeshFilter.h>

namespace rstk {
template<typename TInputMesh, typename TOutputMesh>
class CopyCastMeshFilter: public itk::MeshToMeshFilter<TInputMesh, TOutputMesh> {
public:
	typedef CopyCastMeshFilter Self;
	typedef itk::MeshToMeshFilter<TInputMesh, TOutputMesh> Superclass;
	typedef itk::SmartPointer<Self> Pointer;

	itkNewMacro(Self)
	;itkTypeMacro(CopyCastMeshFilter, itk::MeshToMeshFilter)
	;

	/** Input types. */
	typedef TInputMesh InputMeshType;
	typedef typename InputMeshType::Pointer InputMeshPointer;
	typedef typename InputMeshType::ConstPointer InputMeshConstPointer;
	typedef typename InputMeshType::CoordRepType InputCoordRepType;
	typedef typename InputMeshType::PointType InputPointType;
	typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
	typedef typename InputMeshType::QEPrimal InputQEPrimal;
	typedef typename InputMeshType::VectorType InputVectorType;

	typedef typename InputMeshType::PointDataContainer InputPointDataContainer;
	typedef typename InputMeshType::CellDataContainer InputCellDataContainer;

	typedef typename InputPointDataContainer::ConstPointer InputPointDataContainerConstPointer;
	typedef typename InputMeshType::PointsContainerConstIterator InputPointsContainerConstIterator;
	typedef typename InputMeshType::PointsContainerConstPointer InputPointsContainerConstPointer;
	typedef typename InputMeshType::CellsContainerConstIterator InputCellsContainerConstIterator;
	typedef typename InputMeshType::CellsContainerConstPointer InputCellsContainerConstPointer;

	typedef typename InputMeshType::EdgeCellType InputEdgeCellType;
	typedef typename InputMeshType::PolygonCellType InputPolygonCellType;
	typedef typename InputMeshType::PointIdList InputPointIdList;
	typedef typename InputMeshType::CellTraits InputCellTraits;
	typedef typename InputCellTraits::PointIdInternalIterator InputPointsIdInternalIterator;

	typedef typename InputQEPrimal::IteratorGeom InputQEIterator;

	/** Output types. */
	typedef TOutputMesh OutputMeshType;
	typedef typename OutputMeshType::Pointer OutputMeshPointer;
	typedef typename OutputMeshType::ConstPointer OutputMeshConstPointer;
	typedef typename OutputMeshType::CoordRepType OutputCoordRepType;
	typedef typename OutputMeshType::PointType OutputPointType;
	typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
	typedef typename OutputMeshType::QEPrimal OutputQEPrimal;
	typedef typename OutputMeshType::VectorType OutputVectorType;
	typedef typename OutputQEPrimal::IteratorGeom OutputQEIterator;
	typedef typename OutputMeshType::PointsContainerIterator OutputPointsContainerIterator;
	typedef typename OutputMeshType::PointsContainerPointer OutputPointsContainerPointer;
	typedef typename OutputMeshType::PointsContainerConstPointer OutputPointsContainerConstPointer;

	typedef typename OutputMeshType::PointDataContainer OutputPointDataContainer;
	typedef typename OutputMeshType::CellDataContainer OutputCellDataContainer;
protected:
	CopyCastMeshFilter() {}
	~CopyCastMeshFilter() {}

	void GenerateData() override {
		const InputMeshType *in = this->GetInput();
		OutputMeshType *     out = this->GetOutput();
		// Copy points
		typedef typename TInputMesh::PointsContainerConstPointer InputPointsContainerConstPointer;
		typedef typename TInputMesh::PointsContainerConstIterator InputPointsContainerConstIterator;

		typedef typename TOutputMesh::PointsContainer OutputPointsContainer;
		typedef typename TOutputMesh::PointsContainerPointer OutputPointsContainerPointer;
		typedef typename TOutputMesh::PointType OutputPointType;

		InputPointsContainerConstPointer inPoints = in->GetPoints();

		if (inPoints) {
			InputPointsContainerConstIterator inIt = inPoints->Begin();
			InputPointsContainerConstIterator inEnd = inPoints->End();

			OutputPointsContainerPointer oPoints = out->GetPoints();
			if (oPoints.IsNull()) {
				oPoints = OutputPointsContainer::New();
				out->SetPoints(oPoints);
			}
			OutputPointType pOut;

			while (inIt != inEnd) {
				pOut.CastFrom(inIt.Value());

				// Fix ITK's LPS default orientation
				pOut[0] = -1.0 * pOut[0];
				pOut[1] = -1.0 * pOut[1];
				oPoints->InsertElement(inIt.Index(), pOut);
				++inIt;
			}
		}

		// Copy Edge Cells
		typedef typename TInputMesh::CellsContainer InputCellsContainer;
		typedef typename InputCellsContainer::ConstPointer InputCellsContainerConstPointer;
		typedef typename InputCellsContainer::ConstIterator InputCellsContainerConstIterator;
		typedef typename TInputMesh::EdgeCellType InputEdgeCellType;

		InputCellsContainerConstPointer inEdgeCells = in->GetEdgeCells();

		if (inEdgeCells) {
			InputCellsContainerConstIterator ecIt = inEdgeCells->Begin();
			InputCellsContainerConstIterator ecEnd = inEdgeCells->End();

			while (ecIt != ecEnd) {
				InputEdgeCellType *pe =
						dynamic_cast<InputEdgeCellType *>(ecIt.Value());
				if (pe) {
					out->AddEdgeWithSecurePointList(
							pe->GetQEGeom()->GetOrigin(),
							pe->GetQEGeom()->GetDestination());
				}
				++ecIt;
			}
		}

		// Copy cells
		typedef typename TInputMesh::CellsContainer InputCellsContainer;
		typedef typename InputCellsContainer::ConstPointer InputCellsContainerConstPointer;
		typedef typename InputCellsContainer::ConstIterator InputCellsContainerConstIterator;
		typedef typename TInputMesh::PolygonCellType InputPolygonCellType;
		typedef typename TInputMesh::PointIdList InputPointIdList;
		typedef typename TInputMesh::CellTraits InputCellTraits;
		typedef typename InputCellTraits::PointIdInternalIterator InputPointsIdInternalIterator;

		out->SetCellsAllocationMethod(
				TOutputMesh::CellsAllocatedDynamicallyCellByCell);

		InputCellsContainerConstPointer inCells = in->GetCells();

		if (inCells) {
			InputCellsContainerConstIterator cIt = inCells->Begin();
			InputCellsContainerConstIterator cEnd = inCells->End();
			while (cIt != cEnd) {
				InputPolygonCellType *pe =
						dynamic_cast<InputPolygonCellType *>(cIt.Value());
				if (pe) {
					InputPointIdList points;
					InputPointsIdInternalIterator pIt =
							pe->InternalPointIdsBegin();
					InputPointsIdInternalIterator pEnd =
							pe->InternalPointIdsEnd();

					while (pIt != pEnd) {
						points.push_back((*pIt));
						++pIt;
					}
					out->AddFaceWithSecurePointList(points, false);
				}
				++cIt;
			}
		}

		typedef typename TOutputMesh::PixelType VertexDataType;

		typedef typename TOutputMesh::PointDataContainer OutputPointDataContainer;
		typedef typename OutputPointDataContainer::Pointer OutputPointDataContainerPointer;
		typedef typename TOutputMesh::CellDataContainer OutputCellDataContainer;
		typedef typename OutputCellDataContainer::Pointer OutputCellDataContainerPointer;

		VertexDataType zero = itk::NumericTraits<VertexDataType>::ZeroValue();

		// Set all elements in pointdatacontainer to zero
		OutputPointDataContainerPointer outputPointData =
				OutputPointDataContainer::New();
		size_t nPoints = this->GetInput()->GetPoints()->Size();
		outputPointData->Reserve(nPoints);
		for (size_t i = 0; i < nPoints; i++) {
			outputPointData->SetElement(i, zero);
		}
		this->GetOutput()->SetPointData(outputPointData);

		size_t nCells = this->GetInput()->GetCells()->Size();
		OutputCellDataContainerPointer outCellData =
				OutputCellDataContainer::New();
		outCellData->Reserve(nCells);
		for (size_t i = 0; i < nPoints; i++) {
			outCellData->SetElement(i, zero);
		}
		this->GetOutput()->SetCellData(outputPointData);
	}
};

} //namespace rstk end

#endif /* SOURCE_DIRECTORY__MODULES_FILTERING_INCLUDE_COPYCASTMESHFILTER_H_ */
