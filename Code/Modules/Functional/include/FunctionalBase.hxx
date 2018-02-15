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

#ifndef FUNCTIONALBASE_HXX_
#define FUNCTIONALBASE_HXX_

#include "FunctionalBase.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <memory>
#include <numeric>
#include <vnl/vnl_random.h>
#include <itkImageAlgorithm.h>
#include <itkOrientImageFilter.h>
#include "InternalOrientationFilter.h"
#include <itkContinuousIndex.h>
#include <itkComposeImageFilter.h>

#include <itkIntensityWindowingImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkMeshFileWriter.h>

#include "ComponentsFileWriter.h"

#define MAX_GRADIENT 20.0
#define MIN_GRADIENT 1.0e-8

namespace rstk {

template< typename TReferenceImageType, typename TCoordRepType >
FunctionalBase<TReferenceImageType, TCoordRepType>
::FunctionalBase():
 m_NumberOfContours(0),
 m_NumberOfRegions(2),
 m_NumberOfVertices(0),
 m_SamplingFactor(2),
 m_DecileThreshold(0.05),
 m_DisplacementsUpdated(true),
 m_EnergyUpdated(false),
 m_RegionsUpdated(false),
 m_ApplySmoothing(false),
 m_UseBackground(false),
 m_Value(0.0),
 m_MaxEnergy(0.0)
 {

	this->m_Threader = itk::MultiThreader::New();
	this->m_NumberOfThreads = this->m_Threader->GetNumberOfThreads();

	this->m_Value = itk::NumericTraits<MeasureType>::infinity();
	this->m_Sigma.Fill(0.0);
	this->m_Interp = InterpolatorType::New();
	this->m_MaskInterp = MaskInterpolatorType::New();

	m_InfoBuffer << "{ \"info\": {";
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::Initialize() {
	this->ParseSettings();

	if( this->m_ApplySmoothing ) {
		SmoothingFilterPointer s = SmoothingFilterType::New();
		s->SetInput( this->m_ReferenceImage );
		s->SetSigmaArray( this->m_Sigma );
		s->Update();
		this->m_ReferenceImage = s->GetOutput();
	}

	if (this->m_BackgroundMask.IsNull()) {
		this->m_UseBackground = false;
		ProbabilityMapPointer p = ProbabilityMapType::New();
		p->SetRegions(this->m_ReferenceImage->GetLargestPossibleRegion());
		p->SetSpacing(this->m_ReferenceImage->GetSpacing());
		p->SetOrigin(this->m_ReferenceImage->GetOrigin());
		p->SetDirection(this->m_ReferenceImage->GetDirection());
		p->Allocate();
		p->FillBuffer(0.0);
		this->m_BackgroundMask = p;
	}

	// Initialize corresponding ROI /////////////////////////////
	// Check that high-res reference sampling grid has been initialized
	if ( this->m_ReferenceSamplingGrid.IsNull() ) {
			this->InitializeSamplingGrid();
	}

	// Initialize interpolators
	this->m_Interp->SetInputImage( this->m_ReferenceImage );
	this->m_MaskInterp->SetInputImage(this->m_BackgroundMask);

	// Compute and set regions in m_ROIs
	this->ComputeCurrentRegions();


	// Compute the outer region in each vertex
	this->InitializeContours();

	this->m_RegionValue.SetSize(this->m_NumberOfRegions);
	this->m_RegionValue.Fill(itk::NumericTraits<MeasureType>::ZeroValue());

	this->m_Model = EnergyModelType::New();
	this->m_Model->SetInput(this->m_ReferenceImage);
	this->m_Model->SetMask(this->m_BackgroundMask);
	this->m_Model->SetPriorsMap(this->m_CurrentMaps);
	if(this->m_UseBackground)
		this->m_Model->SetNumberOfSpecialRegions(2);
	this->m_Model->Update();

	this->m_EnergyCalculator = EnergyFilter::New();
	this->m_EnergyCalculator->SetInput(this->m_ReferenceImage);
	this->m_EnergyCalculator->SetPriorsMap(this->m_CurrentMaps);
	this->m_EnergyCalculator->SetMask(this->m_BackgroundMask);
	this->m_EnergyCalculator->SetModel(this->m_Model);
	this->m_EnergyCalculator->Update();

	if( this->m_Priors.size() == this->m_Target.size() ) {
		const MeasureArray finals = this->GetFinalEnergy();
		MeasureType total = 0.0;
		for (size_t i = 0; i < finals.Size(); i++)
			total+=finals[i];

		std::cout << "Energy_total= " << total << " || regions=" << finals << "." << std::endl;
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::LoadShapePriors( std::vector< std::string > movingSurfaceNames ) {
	for( size_t i = 0; i < movingSurfaceNames.size(); i++) {
		typename PriorReader::Pointer polyDataReader = PriorReader::New();
		polyDataReader->SetFileName( movingSurfaceNames[i] );
		polyDataReader->Update();
		this->AddShapePrior( polyDataReader->GetOutput() );
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
size_t
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddShapePrior( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::ScalarContourType* prior ) {
	this->m_Offsets.push_back( this->m_NumberOfVertices );
	this->m_Priors.push_back( prior );

	Scalar2VectorCopyPointer copy = Scalar2VectorCopyType::New();
	copy->SetInput( prior );
	copy->Update();
	this->m_CurrentContours.push_back(copy->GetOutput());

	// Increase number of off-grid nodes to set into the sparse-dense interpolator
	this->m_NumberOfVertices+= prior->GetNumberOfPoints();
	this->m_NumberOfContours++;
	this->m_NumberOfRegions++;

	return this->m_NumberOfContours-1;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeDerivative(PointValueType* grad, ScalesType scales) {
	// Update contours and initialize sizes
	this->UpdateContour();

	size_t nvertices = this->m_ValidVertices.size();

	PointValuesVector gradients;
	gradients.resize(nvertices);
	std::fill(gradients.begin(), gradients.end(), -1.0);

	struct ParallelGradientStruct str;
	str.selfptr = this;
	str.total = nvertices;
	str.gradients = &gradients;


	std::vector<PointDataContainerPointer> normals;
	for (size_t i = 0; i < this->m_NumberOfContours; i++ ) {
		NormalFilterPointer nfilter = NormalFilterType::New();
		nfilter->SetWeight(NormalFilterType::AREA);
		nfilter->SetInput(this->m_CurrentContours[i]);
		nfilter->Update();
		normals.push_back(nfilter->GetOutput()->GetPointData() );

		str.areas.push_back(nfilter->GetVertexAreaContainer());
		str.points.push_back(this->m_CurrentContours[i]->GetPoints());
		str.totalAreas.push_back(nfilter->GetTotalArea());
	}

	// Start multithreading engine
	this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
	this->GetMultiThreader()->SetSingleMethod(this->ThreadedDerivativeCallback, &str);
	this->GetMultiThreader()->SingleMethodExecute();


	PointValuesVector sample(gradients);
	std::sort(sample.begin(), sample.end());

	this->m_GradientStatistics[0] = sample.front();
	this->m_GradientStatistics[1] = sample[int(0.05 * (sample.size()-1))];
	this->m_GradientStatistics[2] = sample[int(0.25 * (sample.size()-1))];
	this->m_GradientStatistics[3] = sample[int(0.50 * (sample.size()-1))];
	this->m_GradientStatistics[4] = sample[int(0.75 * (sample.size()-1))];
	this->m_GradientStatistics[5] = sample[int(0.95 * (sample.size()-1))];
	this->m_GradientStatistics[6] = sample.back();

	VectorType ni, v;
	PointValueType g;
	PointIdentifier pid, cpid;     // id of vertex in its contour
	ROIPixelType icid = 0;
	for(size_t vvid = 0; vvid < nvertices; vvid++ ) {
		pid = this->m_ValidVertices[vvid];
		icid = this->m_InnerRegion[vvid];
		cpid = pid - this->m_Offsets[icid];
		ni = normals[icid]->ElementAt(cpid);
		g = gradients[vvid];
		if ( g > this->m_GradientStatistics[5] ) g = this->m_GradientStatistics[5];
		if ( g < this->m_GradientStatistics[1] ) g = this->m_GradientStatistics[1];

		v.Fill(0.0);
		for( size_t i = 0; i < Dimension; i++ ) {
			if( scales[i] > 1.0e-8 ) {
				v[i] = scales[i] * g * ni[i];
				grad[vvid + i * nvertices] = static_cast<float>(v[i]);
			} else {
				grad[vvid + i * nvertices] = 0.0;
			}
		}

		this->m_CurrentContours[icid]->GetPointData()->SetElement( cpid, v );
	}
}

template< typename TReferenceImageType, typename TCoordRepType >
ITK_THREAD_RETURN_TYPE
FunctionalBase<TReferenceImageType, TCoordRepType>
::ThreadedDerivativeCallback(void *arg) {
	itk::ThreadIdType total, threadId, threadCount;
	threadId = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
	threadCount = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

	ParallelGradientStruct* str = (ParallelGradientStruct *)( ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

	// Compute indices corresponding to each segment
	size_t nvertices = str->total;
	size_t ssize = ceil(1.0 * nvertices / threadCount);
	size_t start = threadId * ssize;
	size_t stop = ( threadId + 1 ) * ssize - 1;

	if (threadId == threadCount - 1)
		stop = nvertices - 1;

	PointValuesVector segment = str->selfptr->ThreadedDerivativeCompute(start, stop, str->points, str->areas, str->totalAreas);

	str->mutex.lock();
	typename PointValuesVector::const_iterator it = segment.begin();
	typename PointValuesVector::const_iterator last = segment.end();
	typename PointValuesVector::iterator dest = str->gradients->begin() + start;

	while(it!=last) {
		*dest = *it;
		++it;
		++dest;
	}
	str->mutex.unlock();

	return ITK_THREAD_RETURN_VALUE;
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointValuesVector
FunctionalBase<TReferenceImageType, TCoordRepType>
::ThreadedDerivativeCompute(size_t start, size_t stop,
		std::vector<PointsContainerPointer> points,
		std::vector<NormalFilterAreasContainer> areas,
		std::vector< double > totalAreas) {
	size_t nvertices = stop - start;
	size_t fullsize = nvertices * Dimension;

	PointIdentifier pid;      // universal id of vertex
	PointIdentifier cpid;     // id of vertex in its contour

	VectorContourPointType ci_prime;
	PointValuesVector sample;
	ROIPixelType ocid;
	ROIPixelType icid = 0;
	double wi = 0.0;

	for(size_t vvid = start; vvid <= stop; vvid++ ) {
		icid = this->m_InnerRegion[vvid];
		ocid = this->m_OuterRegion[vvid];
		pid = this->m_ValidVertices[vvid];
		cpid = pid - this->m_Offsets[icid];
		ci_prime = points[icid]->ElementAt(cpid); // Get c'_i
		wi = areas[icid][cpid] / totalAreas[icid];
		sample.push_back(this->EvaluateGradient( ci_prime, ocid, icid )  * wi);
	}

	return sample;
}


template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::UpdateContour() {
	if( this->m_DisplacementsUpdated ) {
		return;
	}

	MeasureType norm;
	ContinuousIndex point_idx;
	size_t changed = 0;
	size_t gpid = 0;
	std::vector< size_t > invalid;

	std::fill(this->m_OffMaskVertices.begin(), this->m_OffMaskVertices.end(), 0);

	for( size_t contid = 0; contid < this->m_NumberOfContours; contid++ ) {
		typename VectorContourType::PointsContainerConstIterator p_it = this->m_Priors[contid]->GetPoints()->Begin();
		typename VectorContourType::PointsContainerConstIterator p_end = this->m_Priors[contid]->GetPoints()->End();
		PointsContainerPointer curPoints = this->m_CurrentContours[contid]->GetPoints();

		VectorContourPointType ci, ci_prime;
		VectorType disp;
		size_t pid;

		// For all the points in the mesh
		while ( p_it != p_end ) {
			ci = p_it.Value();
			pid = p_it.Index();
			disp = this->m_CurrentDisplacements->GetElement( gpid ); // Get the interpolated value of the field in the point
			norm = disp.GetNorm();

			if( norm > 1.0e-8 ) {
				ci_prime = ci + disp; // Add displacement vector to the point
				if( ! this->CheckExtent(ci_prime,point_idx) ) {
					invalid.push_back( gpid );
					this->InvokeEvent( WarningEvent() );
				}
				curPoints->SetElement( pid, ci_prime );
				changed++;
			}

			if ( (1.0 - this->m_MaskInterp->Evaluate(ci_prime)) < 1.0e-5 ) {
				this->m_OffMaskVertices[contid]++;
			}

			++p_it;
			gpid++;
		}
	}

	if ( invalid.size() > 0 ) {
		itkWarningMacro(<< "a total of " << invalid.size() << " mesh nodes were to be moved off the image domain." );
	}

	this->m_DisplacementsUpdated = true;
	this->m_RegionsUpdated = (changed==0);
	this->m_EnergyUpdated = (changed==0);

	this->ComputeCurrentRegions();
}

template< typename TReferenceImageType, typename TCoordRepType >
typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetValue() {
	if ( !this->m_EnergyUpdated ) {
		this->m_EnergyCalculator->SetPriorsMap(this->m_CurrentMaps);
		this->m_EnergyCalculator->Update();
		this->m_RegionValue = this->m_EnergyCalculator->GetEnergies();

		this->m_Value = 0.0;
		for( size_t roi = 0; roi < this->m_RegionValue.Size(); roi++ ) {
			this->m_Value+= this->m_RegionValue[roi];
		}
		this->m_EnergyUpdated = true;
	}
	return this->m_Value;
}

template< typename TReferenceImageType, typename TCoordRepType >
inline bool
FunctionalBase<TReferenceImageType, TCoordRepType>
::CheckExtent( typename FunctionalBase<TReferenceImageType, TCoordRepType>::VectorContourPointType& p, typename FunctionalBase<TReferenceImageType, TCoordRepType>::ContinuousIndex& idx) const {
	ReferencePointType ref;
	ref.CastFrom ( p );
	bool isInside = this->m_ReferenceImage->TransformPhysicalPointToContinuousIndex( ref , idx );

	if(!isInside) {
		for ( size_t i = 0; i<Dimension; i++) {
			if ( idx[i] < 0.0 ) {
				p.SetElement(i, this->m_FirstPixelCenter[i] );
			}
			else if ( idx[i] > (this->m_ReferenceSize[i] -1) ) {
				p.SetElement(i, this->m_LastPixelCenter[i] );
			}
		}
	}

	return isInside;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeSamplingGrid() {
	typename FieldType::SizeType exp_size;

	for (size_t i = 0; i<Dimension; i++ ){
		exp_size[i] = (unsigned int) (this->m_ReferenceSize[i]*this->m_SamplingFactor);
	}

	PointType firstPixelCenter;
	VectorType step;
	typename FieldType::SpacingType spacing;

	for (size_t i = 0; i<Dimension; i++ ){
		step[i] = (this->m_End[i] - this->m_Origin[i]) / (1.0*exp_size[i]);
		spacing[i]= fabs( step[i] );
		firstPixelCenter[i] = this->m_Origin[i] + 0.5 * step[i];
	}

	this->m_ReferenceSamplingGrid = FieldType::New();
	this->m_ReferenceSamplingGrid->SetOrigin( firstPixelCenter );
	this->m_ReferenceSamplingGrid->SetDirection( this->m_Direction );
	this->m_ReferenceSamplingGrid->SetRegions( exp_size );
	this->m_ReferenceSamplingGrid->SetSpacing( spacing );
	this->m_ReferenceSamplingGrid->Allocate();
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputeCurrentRegions() {
	BinarizeMeshFilterPointer newp = BinarizeMeshFilterType::New();
	newp->SetInputs( this->m_CurrentContours );
	newp->SetOutputReference( this->m_ReferenceSamplingGrid );
	newp->Update();
	this->m_CurrentRegions = newp->GetOutputSegmentation();

	DownsamplePointer p = DownsampleFilter::New();
	p->SetInput(newp->GetOutput());
	p->SetOutputParametersFromImage( this->m_ReferenceImage );
	p->SetMaskImage(m_BackgroundMask);
	p->Update();

	this->m_CurrentMaps = p->GetOutput();
	this->m_RegionsUpdated = true;
}

template< typename TReferenceImageType, typename TCoordRepType >
const typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureArray
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetFinalEnergy() const {
	ScalarContourList groundtruth;

	for(size_t idx = 0; idx < this->m_Target.size(); idx++) {
		ScalarContourCopyPointer copy = ScalarContourCopyType::New();
		copy->SetInput( this->m_Target[idx] );
		copy->Update();
		groundtruth.push_back(copy->GetOutput());
	}

	typedef MultilabelBinarizeMeshFilter< ScalarContourType > ScalarBinarizeMeshFilterType;
	typename ScalarBinarizeMeshFilterType::Pointer newp = ScalarBinarizeMeshFilterType::New();
	newp->SetInputs( groundtruth );
	newp->SetOutputReference( this->m_ReferenceSamplingGrid );
	newp->Update();

	DownsamplePointer p = DownsampleFilter::New();
	p->SetInput(newp->GetOutput());
	p->SetOutputParametersFromImage( this->m_ReferenceImage );
	p->SetMaskImage(m_BackgroundMask);
	p->Update();

	EnergyModelPointer m = EnergyModelType::New();
	m->SetInput(this->m_ReferenceImage);
	m->SetMask(this->m_BackgroundMask);
	m->SetPriorsMap(p->GetOutput());
	if(this->m_UseBackground)
		m->SetNumberOfSpecialRegions(2);
	m->Update();

	EnergyFilterPointer calc = EnergyFilter::New();
	calc->SetInput(this->m_ReferenceImage);
	calc->SetPriorsMap(p->GetOutput());
	calc->SetMask(this->m_BackgroundMask);
	calc->SetModel(m);
	calc->Update();

	const MeasureArray finals = calc->GetEnergies();
	MeasureType total = 0.0;
	for (size_t i = 0; i < finals.Size(); i++)
		total+=finals[i];

	this->m_InfoBuffer << " \"total_energy\": " << total << ", \"roi_energy\": " << finals << ", ";
	this->m_InfoBuffer << " \"model\" : " << m->PrintFormattedDescriptors();
	return finals;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::InitializeContours() {
	this->m_OffMaskVertices.resize(this->m_NumberOfContours);
	std::fill(this->m_OffMaskVertices.begin(), this->m_OffMaskVertices.end(), 0);

	VectorType zerov; zerov.Fill(0.0);
	this->m_CurrentDisplacements = PointDataContainer::New();
	this->m_CurrentDisplacements->Reserve( this->m_NumberOfVertices );

    // Set up ROI interpolator
    typename ROIInterpolatorType::Pointer interp = ROIInterpolatorType::New();
    interp->SetInputImage( this->m_CurrentRegions );

    PointIdentifier tpid = 0;

    // Set up outer regions
    for ( size_t contid = 0; contid < this->m_NumberOfContours; contid ++) {
    	NormalFilterPointer nfilter = NormalFilterType::New();
    	nfilter->SetWeight(NormalFilterType::AREA);
    	nfilter->SetInput(this->m_CurrentContours[contid]);
    	nfilter->Update();
    	PointDataContainer* normals = nfilter->GetOutput()->GetPointData();
    	PointsConstIterator c_it  = this->m_CurrentContours[contid]->GetPoints()->Begin();
    	PointsConstIterator c_end = this->m_CurrentContours[contid]->GetPoints()->End();

    	PointType ci;
    	VectorType v;
    	VectorType ni;
    	size_t pid;
    	float step = 1;
    	while( c_it != c_end ) {

    		pid = c_it.Index();
    		ci = c_it.Value();

    		this->m_Vertices.push_back( ci );
    		this->m_CurrentDisplacements->SetElement(tpid, zerov);

    		if ( (1.0 - this->m_MaskInterp->Evaluate(ci)) < 1.0e-5 ) {
    			this->m_OffMaskVertices[contid]++;
    		}

    		ni = normals->GetElement( pid );
    		ROIPixelType inner = interp->Evaluate( ci + ni );
    		ROIPixelType outer = interp->Evaluate( ci - ni );

    		if (inner == contid) {
    			while (outer==inner) {
    				outer = interp->Evaluate( ci - (ni * step * 0.1) );
    				if (step == 9)
    					break;
    				step++;
    			}

    			if(outer!=inner) {
    				this->m_ValidVertices.push_back(tpid);
    				this->m_OuterRegion.push_back(outer);
    				this->m_InnerRegion.push_back(inner);
    			}
    		}

    		++c_it;
    		tpid++;
    	}
    }
    std::cout << "Valid vertices: " << this->m_ValidVertices.size() << " of " << this->m_Vertices.size() << "." << std::endl;
}

template< typename TReferenceImageType, typename TCoordRepType >
double
FunctionalBase<TReferenceImageType, TCoordRepType>
::ComputePointArea(const PointIdentifier & iId, VectorContourType *mesh ) {
	QEType* edge = mesh->FindEdge( iId );
	QEType* temp = edge;
	CellIdentifier cell_id(0);
	double totalArea = 0.0;
	VectorContourPointType pt[3];
	typedef typename PolygonType::PointIdIterator PolygonPointIterator;

	do {
		cell_id = temp->GetLeft();

		if ( cell_id != VectorContourType::m_NoFace ) {
			PolygonType *poly = dynamic_cast< PolygonType * >(
					mesh->GetCells()->GetElement(cell_id) );
			PolygonPointIterator pit = poly->PointIdsBegin();

			for(size_t k = 0; pit!= poly->PointIdsEnd(); ++pit, k++ ) {
				pt[k] = mesh->GetPoint( *pit );
			}

			totalArea += TriangleType::ComputeArea(pt[0], pt[1], pt[2]);
		}

		temp = temp->GetOnext();
	} while ( temp != edge );
	return fabs(totalArea * 0.33);
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::AddOptions( SettingsDesc& opts ) {
	opts.add_options()
			("smoothing", bpo::value< float > (), "apply isotropic smoothing filter on target image, with kernel sigma=S mm.")
			("smooth-auto", bpo::bool_switch(), "apply isotropic smoothing filter on target image, with automatic computation of kernel sigma.")
			("uniform-bg-membership", bpo::bool_switch(), "consider last ROI as background and do not compute descriptors.")
			("decile-threshold,d", bpo::value< float > (), "set (decile) threshold to consider a computed gradient as outlier (ranges 0.0-0.5)");
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::ParseSettings() {
	if( this->m_Settings.count( "smoothing" ) ) {
		this->m_ApplySmoothing = true;
		bpo::variable_value v = this->m_Settings["smoothing"];
		this->m_Sigma.Fill( v.as<float> () );
	}

	if( this->m_Settings.count( "smooth-auto" ) ) {
		bpo::variable_value v = this->m_Settings["smooth-auto"];
		if ( v.as<bool>() ) {
			this->m_ApplySmoothing= v.as<bool> ();
			this->m_Sigma.Fill( 0.0 );
		}
	}

	if( this->m_Settings.count( "uniform-bg-membership" ) ) {
		bpo::variable_value v = this->m_Settings["uniform-bg-membership"];
		if ( v.as<bool>() ) {
			this->SetUseBackground(true);
		}
	}

	if( this->m_Settings.count( "decile-threshold") ) {
		bpo::variable_value v = this->m_Settings["decile-threshold"];
		this->SetDecileThreshold( v.as<float> () );
	}
	this->Modified();
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::SetCurrentDisplacements( const VNLVectorContainer& vals ) {
	if ( vals.Size() != Dimension ) {
		itkExceptionMacro(<< "vals contains a wrong number of columns");
	}
	size_t npoints = vals[0].size();

	if ( npoints != this->m_NumberOfVertices ) {
		itkExceptionMacro( << "vals contains a wrong number of vectors");
	}

	VectorType new_ci, old_ci;
	MeasureType norm;
	size_t modified = 0;
	for( size_t id = 0; id<npoints; id++ ) {
		new_ci.Fill( 0.0 );
		norm = 0.0;

		for( size_t d = 0; d < Dimension; d++) {
			new_ci[d] = vals[d][id];
		}

		old_ci = this->m_CurrentDisplacements->GetElement( id );
		norm = (new_ci-old_ci).GetNorm();

		if ( norm > 1.0e-8 ) {
			modified++;
			this->m_CurrentDisplacements->SetElement( id, new_ci );
		}
	}

	this->m_DisplacementsUpdated = (modified==0);

}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::LoadReferenceImage ( const std::vector<std::string> fixedImageNames ) {
	typedef itk::ComposeImageFilter< ChannelType, ReferenceImageType >     InputToVectorFilterType;
	typename InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();
	for (size_t i = 0; i < fixedImageNames.size(); i++ ) {
		typename ChannelReader::Pointer r = ChannelReader::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();
		comb->SetInput(i,r->GetOutput());
	}

	comb->Update();
	this->SetReferenceImage(comb->GetOutput());

	// Fix ITK's LPS default orientation
	DirectionType orient; orient.SetIdentity();
	orient[0][0] = -1;
	orient[1][1] = -1;
	this->m_ReferenceImage->SetDirection(orient * this->m_ReferenceImage->GetDirection());
	this->m_ReferenceImage->SetOrigin(orient * this->m_ReferenceImage->GetOrigin());

	// Cache image properties
	this->m_FirstPixelCenter = this->m_ReferenceImage->GetOrigin();
	this->m_Direction = this->m_ReferenceImage->GetDirection();
	this->m_ReferenceSize = this->m_ReferenceImage->GetLargestPossibleRegion().GetSize();
	this->m_ReferenceSpacing = this->m_ReferenceImage->GetSpacing();

	ContinuousIndex tmp_idx;
	tmp_idx.Fill( -0.5 );
	this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( tmp_idx, this->m_Origin );

	for ( size_t dim = 0; dim<FieldType::ImageDimension; dim++)  tmp_idx[dim]= this->m_ReferenceSize[dim]-1.0;
	this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( tmp_idx, this->m_LastPixelCenter );

	for ( size_t dim = 0; dim<FieldType::ImageDimension; dim++)  tmp_idx[dim]= this->m_ReferenceSize[dim]- 0.5;
	this->m_ReferenceImage->TransformContinuousIndexToPhysicalPoint( tmp_idx, this->m_End );

//	std::cout << this->m_FirstPixelCenter << std::endl;
//	std::cout << this->m_LastPixelCenter << std::endl;
//	std::cout << this->m_Direction << std::endl;
}

template< typename TReferenceImageType, typename TCoordRepType >
void
FunctionalBase<TReferenceImageType, TCoordRepType>
::SetBackgroundMask (const ProbabilityMapType * _arg) {
	if ( this->GetDebug() && ::itk::Object::GetGlobalWarningDisplay() ) {
		std::ostringstream itkmsg;
		itkmsg << "Debug: In " "/home/oesteban/workspace/RegSeg/Code/Modules/Functional/include/FunctionalBase.h" ", line " << 333 << "\n" \
				<< this->GetNameOfClass() << " (" << this << "): " "setting " << "BackgroundMask" " to " << _arg       \
				<< "\n\n";
		::itk::OutputWindowDisplayDebugText( itkmsg.str().c_str() );
	}
	if ( this->m_BackgroundMask != _arg ) {
		typedef itk::OrientImageFilter< ProbabilityMapType, ProbabilityMapType >  Orienter;
		typename Orienter::Pointer orient = Orienter::New();
		orient->UseImageDirectionOn();
		orient->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
		orient->SetInput(_arg);
		orient->Update();
		ProbabilityMapPointer msk = orient->GetOutput();
		ReferenceSizeType size = msk->GetLargestPossibleRegion().GetSize();

		DirectionType idmat; idmat.SetIdentity();
		DirectionType itk; itk.SetIdentity();
		itk(0,0) = -1.0; itk(1,1) = -1.0;

		PointType neworig = itk * msk->GetOrigin();
		msk->SetDirection(idmat);
		msk->SetOrigin(neworig);

		typedef itk::IntensityWindowingImageFilter< ProbabilityMapType, ProbabilityMapType > WindowFilter;
		typename WindowFilter::Pointer mskwindow = WindowFilter::New();
		mskwindow->SetInput(msk);
		mskwindow->SetOutputMinimum(0.0);
		mskwindow->SetOutputMaximum(1.0);
		mskwindow->SetWindowMinimum(0.0);
		mskwindow->SetWindowMaximum(1.0);
		mskwindow->Update();

		typedef itk::InvertIntensityImageFilter<ProbabilityMapType> InvertMaskFilterType;
		typename InvertMaskFilterType::Pointer msk_inv = InvertMaskFilterType::New();
		msk_inv->SetInput(mskwindow->GetOutput());
		msk_inv->SetMaximum(1.0);
		msk_inv->Update();

		this->m_BackgroundMask = msk_inv->GetOutput();

		if ( (m_BackgroundMask->GetLargestPossibleRegion() != this->m_ReferenceImage->GetLargestPossibleRegion()) ||
				(m_BackgroundMask->GetOrigin() != this->m_ReferenceImage->GetOrigin())) {
			ProbmapResamplePointer res = ProbmapResampleType::New();
			res->SetInput(msk_inv->GetOutput());
			res->SetSize(this->m_ReferenceSize);
			res->SetOutputSpacing(this->m_ReferenceSpacing);
			res->SetOutputOrigin(this->m_FirstPixelCenter);
			res->SetOutputDirection(this->m_Direction);
			res->SetDefaultPixelValue( 0.0 );
			res->Update();

			this->m_BackgroundMask = res->GetOutput();
		}

		this->Modified();
	}
}


template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::EvaluateGradient( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point,
		size_t outer_roi, size_t inner_roi ) const {
	if ( outer_roi == inner_roi ) {
		return 0.0;
	}
	ReferencePixelType value = this->m_Interp->Evaluate( point );
	MeasureType gin  = this->m_Model->Evaluate( value, inner_roi );
	MeasureType gout = this->m_Model->Evaluate( value, outer_roi );

	MeasureType grad = gin - gout;
	grad = (fabs(grad)>MIN_GRADIENT)?grad:0.0;
	return grad;
}


template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetEnergyAtPoint( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point, size_t roi ) const {
	ReferencePixelType value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	if (isOutside > 1.0)
		isOutside = 1.0;
	else if (isOutside < 1.0e-3)
		isOutside = 0.0;
	float factor = 1.0 * (this->m_NumberOfRegions - roi) / this->m_NumberOfRegions;
	return this->GetEnergyOfSample( value, roi ) + isOutside * factor * this->m_MaxEnergy;
}

template <typename TReferenceImageType, typename TCoordRepType>
inline typename FunctionalBase<TReferenceImageType, TCoordRepType>::MeasureType
FunctionalBase<TReferenceImageType, TCoordRepType>
::GetEnergyAtPoint( const typename FunctionalBase<TReferenceImageType, TCoordRepType>::PointType & point, size_t roi,
		            typename FunctionalBase<TReferenceImageType, TCoordRepType>::ReferencePixelType & value) const {
	value = this->m_Interp->Evaluate( point );
	float isOutside = this->m_MaskInterp->Evaluate( point );
	if (isOutside > 1.0)
		isOutside = 1.0;
	else if (isOutside < 1.0e-3)
		isOutside = 0.0;
	float factor = 1.0 * (this->m_NumberOfRegions - roi) / this->m_NumberOfRegions;
	return this->GetEnergyOfSample( value, roi ) + isOutside * factor * this->m_MaxEnergy;
}

}

#endif /* FUNCTIONALBASE_HXX_ */
