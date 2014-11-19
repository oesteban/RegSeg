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


/**
 * This example generates a sphere, cuts it with a plane and, therefore, generates a circlular contour (vtkPolyData).
 * Subsequently a binary image representation (vtkImageData) is extracted from it. Internally vtkPolyDataToImageStencil and
 * vtkLinearExtrusionFilter are utilized. Both the circular poly data (circle.vtp) and the resultant image (labelImage.mhd)
 * are saved to disk.
 */
int main(int argc, char *argv[]) {
	std::string image;
	std::vector<std::string> surfaces, rsurfaces;
	int slice_num, axis, nimages;
	std::string axisname = "axial";

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
		("help,h", "show help message")
		("input-image,i", bpo::value<std::string>(&image)->required(), "reference image file")
		("surfaces,S", bpo::value<std::vector<std::string>	>(&surfaces)->multitoken(), "surfaces (vtk files)")
		("ref-surfaces,R", bpo::value<std::vector<std::string>	>(&rsurfaces)->multitoken(), "reference surfaces (vtk files)")
		("num-slices,n", bpo::value < int >(&nimages)->default_value(14), "number of slices")
		("slice,s", bpo::value < int >(&slice_num)->default_value(-1), "slice number")
		("axis,a", bpo::value < int >(&axis)->default_value(2), "view");

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


	ImageType::PointType c;
	ImageType::SizeType s = im->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType idx;

	for( size_t i = 0; i < DIMENSION; i++) {
		idx[i] = int(0.5*(s[i]-1));
	}

	if (axis==1) axisname = "coronal";
	if (axis==0) axisname = "sagittal";

	size_t nsurf = surfaces.size();
	size_t nrsurf = rsurfaces.size();


	std::vector< vtkSmartPointer<vtkPolyData> > vp;
	for(size_t surf = 0; surf < surfaces.size(); surf++) {
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(surfaces[surf].c_str());
		reader->Update();
		vp.push_back(reader->GetOutput());
	}

	std::vector< vtkSmartPointer<vtkPolyData> > vpr;
	for(size_t surf = 0; surf < rsurfaces.size(); surf++) {
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(rsurfaces[surf].c_str());
		reader->Update();
		vpr.push_back(reader->GetOutput());
	}

	double totalims = 1.0 * (nimages+1);
	double normal[DIMENSION];
	for (size_t i = 0; i < DIMENSION; i++) {
		normal[i] = (i==axis)?1.0:0.0;
	}

	// Image slice: http://www.cmake.org/Wiki/VTK/Examples/Cxx/Images/ImageSliceMapper
	vtkSmartPointer<vtkImageSliceMapper> imageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
	imageSliceMapper->SliceFacesCameraOn();
	imageSliceMapper->SetOrientation(axis);
	imageSliceMapper->SetInputData(connector->GetOutput());

	for (size_t sl = 0; sl < nimages; sl++) {
		double factor = ((sl+1) / totalims);
		idx[axis] = int( (s[axis]-1)* factor);
		im->TransformIndexToPhysicalPoint(idx, c);

		imageSliceMapper->SetSliceNumber(idx[axis]);

		// Setup renderers
		vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

		for(size_t surf = 0; surf < rsurfaces.size(); surf++) {
			vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
			cutter->SetInputData(vpr[surf]);
			vtkSmartPointer<vtkPlane> cutPlane = vtkSmartPointer<vtkPlane>::New();
			cutPlane->SetOrigin(c[0], c[1], c[2]);
			cutPlane->SetNormal(normal);
			cutter->SetCutFunction(cutPlane);
			vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
			stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
			stripper->Update();

			vtkSmartPointer<vtkPolyDataMapper> contourmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			contourmapper->SetInputConnection(stripper->GetOutputPort());
			vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
			outlineActor->SetMapper(contourmapper);
			outlineActor->GetProperty()->SetColor(1.0,1.0,0.2);
			outlineActor->GetProperty()->SetLineWidth(2);
			renderer->AddActor(outlineActor);
		}

		for(size_t surf = 0; surf < surfaces.size(); surf++) {
			vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
			cutter->SetInputData(vp[surf]);
			vtkSmartPointer<vtkPlane> cutPlane = vtkSmartPointer<vtkPlane>::New();
			cutPlane->SetOrigin(c[0], c[1], c[2]);
			cutPlane->SetNormal(normal);
			cutter->SetCutFunction(cutPlane);
			vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
			stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
			stripper->Update();

			vtkSmartPointer<vtkPolyDataMapper> contourmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			contourmapper->SetInputConnection(stripper->GetOutputPort());
			vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
			outlineActor->SetMapper(contourmapper);
			outlineActor->GetProperty()->SetColor(0.5,0.6,1.0);
			outlineActor->GetProperty()->SetLineWidth(2);
			renderer->AddActor(outlineActor);
		}

		vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
		imageSlice->GetProperty()->SetInterpolationTypeToNearest();
		imageSlice->SetMapper(imageSliceMapper);
		// imageSlice->SetDisplayLocationToBackground();
		renderer->AddViewProp(imageSlice);

		// Setup render window
		int viewExtent[2];
		switch(axis) {
		case 0:
			viewExtent[0] = s[1]*5;
			viewExtent[1] = s[2]*5;
			break;
		case 1:
			viewExtent[0] = s[0]*5;
			viewExtent[1] = s[2]*5;
			break;
		default:
			viewExtent[0] = s[0]*5;
			viewExtent[1] = s[1]*5;
			break;
		}

		std::stringstream slicenumtext;
		slicenumtext << idx[axis];
		vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
		textActor->GetTextProperty()->SetFontSize(40);
		textActor->SetPosition( int(0.25*viewExtent[0]), int(0.25*viewExtent[1]));
		textActor->SetInput( slicenumtext.str().c_str() );
		renderer->AddActor2D( textActor );

//		vtkSmartPointer<vtkCornerAnnotation> cornerAnnotation = vtkSmartPointer<vtkCornerAnnotation>::New();
//		cornerAnnotation->SetLinearFontScaleFactor(2);
//		cornerAnnotation->SetNonlinearFontScaleFactor(1);
//		cornerAnnotation->SetMaximumFontSize(30);
//		cornerAnnotation->SetText( 0, slicenumtext.str().c_str() );
//		renderer->AddViewProp(cornerAnnotation);

		vtkSmartPointer<vtkCamera> cam = renderer->GetActiveCamera();
		cam->SetPosition(normal);
		renderer->ResetCamera();

		vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->SetSize(viewExtent);
		renderWindow->AddRenderer(renderer);

		// Setup render window interactor
	//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//	vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
	//	renderWindowInteractor->SetInteractorStyle(style);
	//	// Render and start interaction
	//	renderWindowInteractor->SetRenderWindow(renderWindow);
	//	renderWindowInteractor->Initialize();
	//
	//	renderWindowInteractor->Start();
	//
		vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
				vtkSmartPointer<vtkWindowToImageFilter>::New();
		windowToImageFilter->SetInput(renderWindow);
		windowToImageFilter->Update();

		vtkSmartPointer<vtkPNGWriter> writer =
				vtkSmartPointer<vtkPNGWriter>::New();

		std::stringstream ss;
		ss << axisname << std::setw(4) << std::setfill('0') << idx[axis] << ".png";
		writer->SetFileName(ss.str().c_str());
		writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		writer->Write();
	}

	return EXIT_SUCCESS;
}
