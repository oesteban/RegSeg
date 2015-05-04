// --------------------------------------------------------------------------------------
// File:          hdistance.cxx
// Date:          Jan 19, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.5.5
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "hdistance.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

int main(int argc,char** argv) {
	std::string refsurf, tstsurf;
	std::string refout = "hdist_ref.vtk";
	std::string tstout = "hdist_tst.vtk";

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
		("help,h", "show help message")
		("reference,r", bpo::value<std::string>(&refsurf)->required(), "reference surface vtk (legacy) file")
		("test,t", bpo::value<std::string>(&tstsurf)->required(), "test surface vtk (legacy) file")
		("reference-output", bpo::value<std::string>(&refout), "reference output surface vtk (legacy) file")
		("test-output", bpo::value<std::string>(&tstout), "test output surface vtk (legacy) file")
		("point-to-cell,C", bpo::bool_switch(), "work in point-to-cell mode (default is point-to-point");

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

	VTK_CREATE(vtkPolyDataReader, rref);
	rref->SetFileName(refsurf.c_str());

	VTK_CREATE(vtkPolyDataReader, rtst);
	rtst->SetFileName(tstsurf.c_str());

	VTK_CREATE(vtkHausdorffDistancePointSetFilter, filter);
	filter->SetInputConnection(rref->GetOutputPort());
	filter->SetInputConnection(1, rtst->GetOutputPort());

	if (vm.count("point-to-cell") && vm["point-to-cell"].as<bool>())
		filter->SetTargetDistanceMethod(1);
	filter->Update();

	VTK_CREATE(vtkPolyDataWriter, polyDataWriter);
	polyDataWriter->SetInputConnection(filter->GetOutputPort(0));
	polyDataWriter->SetFileName(refout.c_str());
	polyDataWriter->Update();

	polyDataWriter->SetInputConnection(filter->GetOutputPort(1));
	polyDataWriter->SetFileName(tstout.c_str());
	polyDataWriter->Update();

	return (EXIT_SUCCESS);
}
