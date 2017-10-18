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
