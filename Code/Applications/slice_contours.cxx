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
	std::vector<std::string> surfaces;
	int slice_num, axis;
	int nimages = 8;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
		("help,h", "show help message")
		("input-image,i", bpo::value<std::string>(&image)->required(), "reference image file")
		("surfaces,S", bpo::value<std::vector<std::string>	>(&surfaces)->multitoken(), "surfaces (vtk files)")
		("slice,s", bpo::value < int >(&slice_num)->default_value(-1), "slice number")
		("plane,p", bpo::value < int >(&axis)->default_value(0), "view");

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

	double totalims = 1.0 * (nimages+1);
	for (size_t sl = 0; sl < nimages; sl++) {
		double factor = ((sl+1) / totalims);
		idx[2] = int( (s[2]-1)* factor);
		std::cout << factor << ", " << idx[2] << std::endl;
		im->TransformIndexToPhysicalPoint(idx, c);

		// Image slice: http://www.cmake.org/Wiki/VTK/Examples/Cxx/Images/ImageSliceMapper
		vtkSmartPointer<vtkImageSliceMapper> imageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
		imageSliceMapper->SetSliceNumber(idx[2]);
		imageSliceMapper->SetInputData(connector->GetOutput());

		// Setup renderers
		vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();


		for(size_t surf = 0; surf < surfaces.size(); surf++) {
			vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(surfaces[surf].c_str());
			reader->Update();
			vtkSmartPointer<vtkPolyData> vp = reader->GetOutput();
			vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
			cutter->SetInputData(vp);
			vtkSmartPointer<vtkPlane> cutPlane = vtkSmartPointer<vtkPlane>::New();
			cutPlane->SetOrigin(c[0], c[1], c[2]);
			cutPlane->SetNormal(0, 0, 1);
			cutter->SetCutFunction(cutPlane);
			vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
			stripper->SetInputConnection(cutter->GetOutputPort()); // valid circle
			stripper->Update();

			vtkSmartPointer<vtkPolyDataMapper> contourmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			contourmapper->SetInputConnection(stripper->GetOutputPort());
			vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
			outlineActor->SetMapper(contourmapper);
			outlineActor->GetProperty()->SetColor(1.0,0,0);
			outlineActor->GetProperty()->SetLineWidth(2);
			renderer->AddActor(outlineActor);
		}

		vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
		imageSlice->GetProperty()->SetInterpolationTypeToNearest();
		imageSlice->SetMapper(imageSliceMapper);
		// imageSlice->SetDisplayLocationToBackground();
		renderer->AddViewProp(imageSlice);
		renderer->ResetCamera();

		// Setup render window
		vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindow->SetSize(s[0]*5, s[1]*5);
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
		ss << "axial" << std::setw(4) << std::setfill('0') << idx[2] << ".png";
		writer->SetFileName(ss.str().c_str());
		writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		writer->Write();
	}

	return EXIT_SUCCESS;
}
