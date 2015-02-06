// --------------------------------------------------------------------------------------
// File:          rbf_approx.cxx
// Date:          Feb 5, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "rbf_approx.h"

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;


int main(int argc, char *argv[]) {
	std::string outPrefix = "";
	std::string maskname, inputname;
	std::vector<size_t> grid_size;
	std::vector<float> grid_spacing;

	bool hasMask = false;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("input,i", bpo::value< std::string >(&inputname)->required(), "input image")
			("mask,m", bpo::value< std::string >(&maskname), "input image mask")
			("grid-size,g", bpo::value< std::vector<size_t> >(&grid_size)->multitoken(), "size of grid of RBF control points")
			("grid-spacing,s", bpo::value< std::vector<float> >(&grid_spacing)->multitoken(), "spacing (mm) between RBF control points");

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

	TransformPointer tf = TransformType::New();

	ReaderPointer r = ReaderType::New();
	r->SetFileName(inputname);
	r->Update();

	ImagePointer im = r->GetOutput();
	tf->SetDomainExtent(im);

	ImagePointer msk;
	if (vm.count("mask")) {
		ReaderPointer rm = ReaderType::New();
		rm->SetFileName(maskname);
		rm->Update();
		msk = rm->GetOutput();
		hasMask = true;
	}

	if (vm.count("grid-spacing")) {
		bpo::variable_value v = vm["grid-spacing"];
		std::vector<float> s = v.as< std::vector<float> > ();
		typename TransformType::SpacingType sp;
		if (s.size() == 1) {
			sp.Fill(s[0]);
		} else if (s.size() == Dimension) {
			for( size_t i = 0; i < Dimension; i++)
				sp[i] = s[i];
		}

		tf->SetControlGridSpacing(sp);
	} else if (vm.count("grid-size")) {
		bpo::variable_value v = vm["grid-size"];
		std::vector<size_t> s = v.as< std::vector<size_t> > ();
		typename TransformType::SizeType size;
		if (s.size() == 1) {
			size.Fill(s[0]);
		} else if (s.size() == Dimension) {
			for( size_t i = 0; i < Dimension; i++)
				size[i] = s[i];
		}

		tf->SetControlGridSize(size);
	} else {
		std::cerr << "Error: at least one of --grid-spacing or --grid-size should be specified" << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	tf->SetOutputReference(im);

	PixelType* m = (hasMask)?msk->GetBufferPointer():NULL;
	PixelType* v = im->GetBufferPointer();
	size_t nPix = im->GetLargestPossibleRegion().GetNumberOfPixels();

	float mval = 1.0;
	float val = 0.0;
	for (size_t i = 0; i < nPix; i++) {
		if(hasMask) {
			mval = *(m + i);
		}

		if (mval > 0.0) {
			val = *(v + i);
			tf->SetPointValue(i, val);
		}
	}

	tf->Initialize();
	tf->ComputeCoefficients();
	tf->Interpolate();

	for (size_t i = 0; i < nPix; i++) {
		*(v + i) = tf->GetPointValue(i)[0];
	}

	WriterPointer w = WriterType::New();
	w->SetFileName("test.nii.gz");
	w->SetInput(im);
	w->Update();

	return EXIT_SUCCESS;
} // main
