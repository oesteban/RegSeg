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
	std::vector<std::string> surfaces, rsurfaces;
	int slice_num, nimages;
	std::string axisname = "axial";
	std::vector<int> axislist;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
		("help,h", "show help message")
		("input-image,i", bpo::value<std::string>(&image)->required(), "reference image file")
		("surfaces,S", bpo::value<std::vector<std::string>	>(&surfaces)->multitoken(), "surfaces (vtk files)")
		("ref-surfaces,R", bpo::value<std::vector<std::string>	>(&rsurfaces)->multitoken(), "reference surfaces (vtk files)")
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

	// Setup render window
	int viewExtent[2] = { 600, 600 };

	// Setup renderers
	VTK_CREATE(vtkRenderer, renderer);

	VTK_CREATE(vtkImageReslice, reslice);
	reslice->SetOutputDimensionality(2);
	reslice->InterpolateOff();
	reslice->SetInputData(connector->GetOutput());
	reslice->SetOutputSpacing( vtkim->GetSpacing() );

	VTK_CREATE(vtkLookupTable, table);
	table->SetRange(0, 255); // image intensity range
	table->SetValueRange(0.0, 1.0); // from black to white
	table->SetSaturationRange(0.0, 0.0); // no color saturation
	table->SetRampToLinear();
	table->Build();

	VTK_CREATE(vtkImageMapToColors, color);
	color->SetLookupTable(table);
	color->SetInputConnection(reslice->GetOutputPort());

	VTK_CREATE(vtkImageActor, imageActor);
	imageActor->GetMapper()->SetInputConnection(color->GetOutputPort());
	renderer->AddActor(imageActor);



	//VTK_CREATE(vtkImageMagnify, imgMagnify);
	//imgMagnify->SetInputData(connector->GetOutput());
	//imgMagnify->SetMagnificationFactors(5, 5, 5);
	// VTK_CREATE(vtkImageViewer2, imageViewer);
	// imageViewer->SetInputData(connector->GetOutput());
	//imageViewer->SetInputConnection(imgMagnify->GetOutputPort());
	// imageViewer->GetImageActor()->InterpolateOff();
	// renderer->AddActor( imageViewer->GetImageActor() );

	VTK_CREATE(vtkTextActor, textActor);
	textActor->GetTextProperty()->SetFontSize(40);
	renderer->AddActor2D( textActor );

	VTK_CREATE(vtkRenderWindow, renderWindow);
	// renderWindow->SetOffScreenRendering(1);
	//renderWindow->SetSize(viewExtent);
	// renderWindow->OffScreenRenderingOn();
	renderWindow->AddRenderer(renderer);



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
		axislist.empty();

		for (size_t aa = 0; aa < 3; aa++) axislist.push_back(aa);
	}

	size_t nsurf = surfaces.size();
	size_t nrsurf = rsurfaces.size();

	std::vector< vtkSmartPointer<vtkPolyData> > vp;
	for(size_t surf = 0; surf < surfaces.size(); surf++) {
		VTK_CREATE(vtkPolyDataReader, reader);
		reader->SetFileName(surfaces[surf].c_str());
		reader->Update();
		vp.push_back(reader->GetOutput());
	}

	std::vector< vtkSmartPointer<vtkPolyData> > vpr;
	for(size_t surf = 0; surf < rsurfaces.size(); surf++) {
		VTK_CREATE(vtkPolyDataReader, reader);
		reader->SetFileName(rsurfaces[surf].c_str());
		reader->Update();
		vpr.push_back(reader->GetOutput());
	}

	double totalims = 1.0 * (nimages+1);

	for(size_t ax = 0; ax < axislist.size(); ax++) {
		int axis = axislist[ax];
		if (axis==1) axisname = "coronal";
		if (axis==0) axisname = "sagittal";

		double normal[DIMENSION];
		for (size_t i = 0; i < DIMENSION; i++) {
			normal[i] = (i==axis)?1.0:0.0;
		}


		// imageViewer->SetSliceOrientation(axis);
		// Image slice: http://www.cmake.org/Wiki/VTK/Examples/Cxx/Images/ImageSliceMapper
		// VTK_CREATE(vtkImageSliceMapper, imageSliceMapper);
		// imageSliceMapper->SliceFacesCameraOn();
		// imageSliceMapper->SetOrientation(axis);
		// imageSliceMapper->SetInputData(connector->GetOutput());

		for (size_t sl = 0; sl < nimages; sl++) {
			VTK_CREATE(vtkMatrix4x4, mat);

			double factor = ((sl+1) / totalims);
			idx[axis] = int( (s[axis]-1)* factor);
			im->TransformIndexToPhysicalPoint(idx, c);

			int slice2DExtent[4];
			int sliceExtent[6];
			vtkim->GetExtent(sliceExtent);
			double sliceOrigin[3];
			vtkim->GetOrigin(sliceOrigin);

			switch(axis) {
			case 0:
				mat->DeepCopy(sagittalElements);
				mat->SetElement(0, 3, c[0]);
			    slice2DExtent[0] = sliceExtent[2];
			    slice2DExtent[1] = sliceExtent[3];
			    slice2DExtent[2] = sliceExtent[4];
			    slice2DExtent[3] = sliceExtent[5];
			    sliceExtent[0] = 0;
			    sliceExtent[1] = 0;
			    sliceOrigin[0] = c[0];
			    break;
				break;
			case 1:
				mat->DeepCopy(coronalElements);
				mat->SetElement(0, 3, c[1]);
			    slice2DExtent[0] = sliceExtent[0];
			    slice2DExtent[1] = sliceExtent[1];
			    slice2DExtent[2] = sliceExtent[4];
			    slice2DExtent[3] = sliceExtent[5];
			    sliceExtent[2] = 0;
			    sliceExtent[3] = 0;
			    sliceOrigin[1] = c[1];
				break;
			default:
				mat->DeepCopy(axialElements);
				mat->SetElement(0, 3, c[2]);
			    slice2DExtent[0] = sliceExtent[0];
			    slice2DExtent[1] = sliceExtent[1];
			    slice2DExtent[2] = sliceExtent[2];
			    slice2DExtent[3] = sliceExtent[3];
			    sliceExtent[4] = 0;
			    sliceExtent[5] = 0;
			    sliceOrigin[2] = c[2];
				break;
			}

			std::cout << "Normal=(" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")";
			std::cout << ", origin=(" << sliceOrigin[0] << ", " << sliceOrigin[1] << ", " << sliceOrigin[2] << ")";
			std::cout <<", extent=[" << sliceExtent[0] << ", " << sliceExtent[1] << ", " << sliceExtent[2] << ", ";
			std::cout <<", " << sliceExtent[3] << ", " << sliceExtent[4] << ", " << sliceExtent[5] << "]." << std::endl;

			reslice->SetResliceAxes(mat);
			reslice->SetOutputExtent(slice2DExtent);
			reslice->Update();

			// imageSliceMapper->SetSliceNumber(idx[axis]);
			// imageViewer->SetSlice(idx[axis]);


//			for(size_t surf = 0; surf < rsurfaces.size(); surf++) {
//				VTK_CREATE(vtkCutter, cutter);
//#if VTK_MAJOR_VERSION <= 5
//				cutter->SetInput(vpr[surf]);
//#else
//				cutter->SetInputData(vpr[surf]);
//#endif
//				VTK_CREATE(vtkPlane, cutPlane);
//				cutPlane->SetOrigin(c[0], c[1], c[2]);
//				cutPlane->SetNormal(normal);
//				cutter->SetCutFunction(cutPlane);
//				VTK_CREATE(vtkStripper, stripper);
//				stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
//				stripper->Update();
//
//				VTK_CREATE(vtkPolyDataMapper, contourmapper);
//				contourmapper->SetInputConnection(stripper->GetOutputPort());
//				VTK_CREATE(vtkActor, outlineActor);
//				outlineActor->SetMapper(contourmapper);
//				outlineActor->GetProperty()->SetColor(1.0,1.0,0.2);
//				outlineActor->GetProperty()->SetLineWidth(2);
//				renderer->AddActor(outlineActor);
//			}
//
//			for(size_t surf = 0; surf < surfaces.size(); surf++) {
//				VTK_CREATE(vtkCutter, cutter);
//#if VTK_MAJOR_VERSION <= 5
//				cutter->SetInput(vp[surf]);
//#else
//				cutter->SetInputData(vp[surf]);
//#endif
//				VTK_CREATE(vtkPlane, cutPlane);
//				cutPlane->SetOrigin(c[0], c[1], c[2]);
//				cutPlane->SetNormal(normal);
//				cutter->SetCutFunction(cutPlane);
//				VTK_CREATE(vtkStripper, stripper);
//				stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
//				stripper->Update();
//				VTK_CREATE(vtkPolyDataMapper, contourmapper);
//				contourmapper->SetInputConnection(stripper->GetOutputPort());
//				VTK_CREATE(vtkActor, outlineActor);
//				outlineActor->SetMapper(contourmapper);
//				outlineActor->GetProperty()->SetColor(0.5,0.6,1.0);
//				outlineActor->GetProperty()->SetLineWidth(2);
//				renderer->AddActor(outlineActor);
//			}

			// VTK_CREATE(vtkImageSlice, imageSlice);
			// imageSlice->GetProperty()->SetInterpolationTypeToNearest();
			// imageSlice->SetMapper(imageSliceMapper);

			// renderer->AddViewProp(imageSlice);
			//renderer->AddRenderer

			std::stringstream slicenumtext;
			slicenumtext << idx[axis];

			textActor->SetPosition( int(0.25*viewExtent[0]), int(0.25*viewExtent[1]));
			textActor->SetInput( slicenumtext.str().c_str() );

//			vtkSmartPointer<vtkCamera> cam = renderer->GetActiveCamera();
//			cam->SetPosition(normal);
//			renderer->ResetCamera();

			VTK_CREATE(vtkWindowToImageFilter, windowToImageFilter);
			windowToImageFilter->SetInput(renderWindow);
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
