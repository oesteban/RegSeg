// --------------------------------------------------------------------------------------
// File:          slice_contours.cxx
// Date:          Nov 18, 2014
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "slice_contours.h"


// see: https://github.com/InsightSoftwareConsortium/LesionSizingToolkit/blob/master/Utilities/Visualization/ViewImageSlicesAndSegmentationContours.cxx
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

/**
 * This example generates a sphere, cuts it with a plane and, therefore, generates a circlular contour (vtkPolyData).
 * Subsequently a binary image representation (vtkImageData) is extracted from it. Internally vtkPolyDataToImageStencil and
 * vtkLinearExtrusionFilter are utilized. Both the circular poly data (circle.vtp) and the resultant image (labelImage.mhd)
 * are saved to disk.
 */
int main(int argc, char *argv[]) {

#if VTK_MAJOR_VERSION <= 5
	// Setup offscreen rendering
	VTK_CREATE(vtkGraphicsFactory, graphics_factory);
	graphics_factory->SetOffScreenOnlyMode( 1);
	graphics_factory->SetUseMesaClasses( 1 );

	VTK_CREATE(vtkImagingFactory, imaging_factory);
	imaging_factory->SetUseMesaClasses( 1 );
#endif

	std::string image;
	std::vector<std::string> surfaces, rsurfaces, osurfaces;
	int slice_num, nimages;
	std::string axisname = "axial";
	std::vector<int> axislist;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
		("help,h", "show help message")
		("input-image,i", bpo::value<std::string>(&image)->required(), "reference image file")
		("surfaces,S", bpo::value<std::vector<std::string>	>(&surfaces)->multitoken(), "surfaces (vtk files)")
		("ref-surfaces,R", bpo::value<std::vector<std::string>	>(&rsurfaces)->multitoken(), "reference surfaces (vtk files)")
		("orig-surfaces,O", bpo::value<std::vector<std::string>	>(&osurfaces)->multitoken(), "original surfaces (vtk files)")
		("num-slices,n", bpo::value < int >(&nimages)->default_value(14), "number of slices")
		("slice,s", bpo::value < int >(&slice_num)->default_value(-1), "slice number")
		("axis,a", bpo::value < std::vector<int> >(&axislist), "axes to be extracted")
		("all-axis,A", bpo::bool_switch(), "export all axes");

	bpo::variables_map vm;

	try {
		bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);

		if (vm.count("help")) {
			std::cout << all_desc << std::endl;
			return 1;
		}

		bpo::notify(vm);
	} catch ( bpo::error& e ) {
		std::cerr << "Error: " << e.what() << std::endl << std::endl;
		std::cerr << all_desc << std::endl;
		return 1;
	}

	ReaderPointer r = ReaderType::New();
	r->SetFileName(image);
	r->Update();

	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(r->GetOutput());
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();

	ImagePointer im = rescaleFilter->GetOutput();
	ImageType::DirectionType dir;
	dir.SetIdentity();
	ImageType::PointType zero; zero.Fill(0.0);
	ImageType::PointType orig = im->GetOrigin();
	orig = im->GetDirection() * (orig - zero);
	im->SetOrigin(orig);
	im->SetDirection(dir);

	ConnectorType::Pointer connector = ConnectorType::New();
	connector->SetInput(im);
	connector->Update();
	vtkSmartPointer<vtkImageData> vtkim = connector->GetOutput();

	double spacing[3];
	vtkim->GetSpacing(spacing);

	// Setup render window
	int viewExtent[2] = { 600, 600 };


	ImageType::PointType c;
	ImageType::SizeType s = im->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType idx;

	for( size_t i = 0; i < DIMENSION; i++) {
		idx[i] = int(0.5*(s[i]-1));
	}

	if (axislist.size() == 0) {
		axislist.push_back(2);
	}

	if (vm.count("all-axis") && vm["all-axis"].as<bool>()) {
		axislist.resize(0);

		for (size_t aa = 0; aa < 3; aa++) axislist.push_back(aa);
	}

	size_t nsurf = surfaces.size();
	size_t nrsurf = rsurfaces.size();
	size_t nosurf = osurfaces.size();

	std::vector< vtkSmartPointer<vtkPolyData> > vp;
	for(size_t surf = 0; surf < surfaces.size(); surf++) {
		VTK_CREATE(vtkPolyDataReader, reader);
		reader->SetFileName(surfaces[surf].c_str());
		reader->Update();
		vp.push_back(reader->GetOutput());
	}

	std::vector< vtkSmartPointer<vtkPolyData> > vpr;
	for(size_t surf = 0; surf < nrsurf; surf++) {
		VTK_CREATE(vtkPolyDataReader, reader);
		reader->SetFileName(rsurfaces[surf].c_str());
		reader->Update();
		vpr.push_back(reader->GetOutput());
	}

	std::vector< vtkSmartPointer<vtkPolyData> > vpo;
	for(size_t surf = 0; surf < nosurf; surf++) {
		VTK_CREATE(vtkPolyDataReader, reader);
		reader->SetFileName(osurfaces[surf].c_str());
		reader->Update();
		vpo.push_back(reader->GetOutput());
	}

	double totalims = 1.0 * (nimages+1);
	int slice2DExtent[4];
	int sliceExtent[6];
	double sliceOrigin[3];
	float height = 0.0;

	for(size_t ax = 0; ax < axislist.size(); ax++) {
		int axis = axislist[ax];
		vtkim->GetExtent(sliceExtent);
		VTK_CREATE(vtkMatrix4x4, mat);
		switch(axis) {
		case 0:
			axisname = "sagittal";
			mat->DeepCopy(sagittalElements);
			mat->SetElement(0, 3, c[0]);
		    height = sliceExtent[1];
		    slice2DExtent[0] = sliceExtent[2];
		    slice2DExtent[1] = sliceExtent[3];
		    slice2DExtent[2] = sliceExtent[4];
		    slice2DExtent[3] = sliceExtent[5];
		    sliceExtent[0] = 0;
		    sliceExtent[1] = 0;
		    break;
		case 1:
			axisname = "coronal";
			mat->DeepCopy(coronalElements);
			mat->SetElement(0, 3, c[1]);
			height = sliceExtent[3];
		    slice2DExtent[0] = sliceExtent[0];
		    slice2DExtent[1] = sliceExtent[1];
		    slice2DExtent[2] = sliceExtent[4];
		    slice2DExtent[3] = sliceExtent[5];
		    sliceExtent[2] = 0;
		    sliceExtent[3] = 0;
			break;
		case 2:
			axisname = "axial";
			mat->DeepCopy(axialElements);
			mat->SetElement(0, 3, c[2]);
			height = sliceExtent[5];
		    slice2DExtent[0] = sliceExtent[0];
		    slice2DExtent[1] = sliceExtent[1];
		    slice2DExtent[2] = sliceExtent[2];
		    slice2DExtent[3] = sliceExtent[3];
		    sliceExtent[4] = 0;
		    sliceExtent[5] = 0;
			break;
		}

		double normal[DIMENSION] = { 0.0, 0.0, 0.0 };
		normal[axis] = 1.0;

		// std::cout << "Axis (" << axisname << ") Normal[" << axis<< "]=[" << normal[0] << ", " << normal[1] << ", " << normal[2] << "]";

		for (size_t sl = 0; sl < nimages; sl++) {
			// Setup renderers
			VTK_CREATE(vtkRenderer, renderer);
			renderer->SetViewport(0, 0, 1, 1);
			renderer->SetBackground(1.0, 1.0, 1.0);
			VTK_CREATE(vtkRenderWindow, renderWindow);
			//renderWindow->SetSize(viewExtent);
			renderWindow->OffScreenRenderingOn();
			renderWindow->Render();
			renderWindow->AddRenderer(renderer);

			double factor = ((sl+1) / totalims);
			idx[axis] = int( (s[axis]-1)* factor);
			im->TransformIndexToPhysicalPoint(idx, c);

			vtkim->GetOrigin(sliceOrigin);
			sliceOrigin[axis] = c[axis];

			VTK_CREATE(vtkPlane, cutPlane);
			cutPlane->SetOrigin(sliceOrigin);
			cutPlane->SetNormal(normal);


			// Setup image slice mapper
			VTK_CREATE(vtkImageResliceMapper, imSliceMap);
			imSliceMap->SetInputData(connector->GetOutput());
			imSliceMap->BorderOn();
			imSliceMap->SetSlicePlane(cutPlane);
			imSliceMap->ResampleToScreenPixelsOff();
			imSliceMap->SliceAtFocalPointOff();
			imSliceMap->UpdateInformation();

			// Image slice: http://www.cmake.org/Wiki/VTK/Examples/Cxx/Images/ImageSliceMapper
			VTK_CREATE(vtkImageSlice, imSlice);
			imSlice->SetMapper(imSliceMap);
			imSlice->GetProperty()->SetInterpolationTypeToNearest();
			renderer->AddViewProp(imSlice);

			for(size_t surf = 0; surf < rsurfaces.size(); surf++) {
				VTK_CREATE(vtkCutter, cutter);
				cutter->SetCutFunction(cutPlane);
#if VTK_MAJOR_VERSION <= 5
				cutter->SetInput(vpr[surf]);
#else
				cutter->SetInputData(vpr[surf]);
#endif

				VTK_CREATE(vtkStripper, stripper);
				stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
				stripper->Update();

				VTK_CREATE(vtkPolyDataMapper, contourmapper);
				contourmapper->SetInputConnection(stripper->GetOutputPort());
				VTK_CREATE(vtkActor, outlineActor);
				outlineActor->SetMapper(contourmapper);
				outlineActor->GetProperty()->SetColor(0.0, 0.804, 0.0);
				outlineActor->GetProperty()->SetLineWidth(4);
				renderer->AddActor(outlineActor);
			}

			for(size_t surf = 0; surf < osurfaces.size(); surf++) {
				VTK_CREATE(vtkCutter, cutter);
				cutter->SetCutFunction(cutPlane);
#if VTK_MAJOR_VERSION <= 5
				cutter->SetInput(vpo[surf]);
#else
				cutter->SetInputData(vpo[surf]);
#endif

				VTK_CREATE(vtkStripper, stripper);
				stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
				stripper->Update();

				VTK_CREATE(vtkPolyDataMapper, contourmapper);
				contourmapper->SetInputConnection(stripper->GetOutputPort());
				VTK_CREATE(vtkActor, outlineActor);
				outlineActor->SetMapper(contourmapper);
				outlineActor->GetProperty()->SetColor(0.118, 0.565, 1.0);
				outlineActor->GetProperty()->SetLineWidth(1.5);
				renderer->AddActor(outlineActor);
			}

			for(size_t surf = 0; surf < surfaces.size(); surf++) {
				VTK_CREATE(vtkCutter, cutter);
	            cutter->SetCutFunction(cutPlane);
#if VTK_MAJOR_VERSION <= 5
				cutter->SetInput(vp[surf]);
#else
				cutter->SetInputData(vp[surf]);
#endif
				VTK_CREATE(vtkStripper, stripper);
				stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
				stripper->Update();
				VTK_CREATE(vtkPolyDataMapper, contourmapper);
				contourmapper->SetInputConnection(stripper->GetOutputPort());
				VTK_CREATE(vtkActor, outlineActor);
				outlineActor->SetMapper(contourmapper);
				outlineActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
				outlineActor->GetProperty()->SetLineWidth(4);
				renderer->AddActor(outlineActor);
			}

			std::stringstream slicenumtext;
			slicenumtext << idx[axis];


			VTK_CREATE(vtkTextActor, textActor);
			textActor->GetTextProperty()->SetFontSize(40);
			renderer->AddActor2D( textActor );
			textActor->SetPosition( 10, 10 );
			textActor->SetInput( slicenumtext.str().c_str() );

//			VTK_CREATE(vtkTextActor, loc1);
//			loc1->GetTextProperty()->SetFontSize(40);
//			renderer->AddActor2D( loc1 );
//			loc1->SetPosition( 280, 10 );
//			loc1->SetInput( std::string("I").c_str() );
//
//			VTK_CREATE(vtkTextActor, loc2);
//			loc2->GetTextProperty()->SetFontSize(40);
//			renderer->AddActor2D( loc2 );
//			loc2->SetPosition( 280, 550 );
//			loc2->SetInput( std::string("S").c_str() );


			vtkSmartPointer<vtkCamera> cam = renderer->GetActiveCamera();
			cam->SetPosition(normal);
			renderer->ResetCamera();
			cam->ParallelProjectionOn();
			cam->SetParallelScale(0.5 * spacing[axis] * height);
			renderWindow->Render();

//			VTK_CREATE(vtkCubeAxesActor2D, axes);
//			axes->SetViewProp(imSlice);
//			axes->SetCamera(cam);

			VTK_CREATE(vtkWindowToImageFilter, windowToImageFilter);
			windowToImageFilter->SetInput(renderWindow);
			windowToImageFilter->SetMagnification(2);
			windowToImageFilter->SetInputBufferTypeToRGBA();
			windowToImageFilter->ReadFrontBufferOff();
			windowToImageFilter->Update();

			VTK_CREATE(vtkPNGWriter, writer);
			writer->SetInputConnection(windowToImageFilter->GetOutputPort());

			std::stringstream ss;
			ss << axisname << std::setw(4) << std::setfill('0') << idx[axis] << ".png";
			writer->SetFileName(ss.str().c_str());
			writer->Write();
		}
	}

	return EXIT_SUCCESS;
}
