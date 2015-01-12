// --------------------------------------------------------------------------------------
// File:          binarize_surfaces.cxx
// Date:          Jan 12, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "binarize_surfaces.h"

int main(int argc, char *argv[]) {
	std::string outPrefix = "binarized";
	std::string refname;
	std::vector< std::string > surfnames;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("surfaces,S", bpo::value < std::vector< std::string > >(&surfnames)->multitoken()->required(), "input surfaces" )
			("reference,R", bpo::value< std::string >(&refname)->required(), "reference image" )
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files");

	bpo::variables_map vm;

	try {
		bpo::store(	bpo::parse_command_line( argc, argv, all_desc ), vm);
		if (vm.count("help") || vm.size() == 0 ) {
			std::cout << all_desc << std::endl;
			return 1;
		}
		bpo::notify(vm);
	} catch ( bpo::error& e ) {
		std::cerr << "Error: " << e.what() << std::endl << std::endl;
		std::cerr << all_desc << std::endl;
		return 1;
	}

	InputMeshContainer priors;

	for( size_t i = 0; i<surfnames.size(); i++) {
		typename ReaderType::Pointer r = ReaderType::New();
		r->SetFileName(surfnames[i]);
		r->Update();
		priors.push_back(r->GetOutput());
	}

	typename ImageReader::Pointer rr = ImageReader::New();
	rr->SetFileName(refname);
	rr->Update();

	// TODO: fix input orientation here.


	BinarizeMeshFilterPointer newp = BinarizeMeshFilterType::New();
	newp->SetInputs( priors );
	newp->SetOutputReference( rr->GetOutput() );
	newp->Update();

	std::stringstream ss;
	ss << outPrefix << "_segs";
	typename SegmentationWriter::Pointer sw = SegmentationWriter::New();
	sw->SetFileName(ss.str().c_str());
	sw->SetInput(newp->GetOutputSegmentation());
	sw->Update();

	DownsamplePointer p = DownsampleFilter::New();
	p->SetInput(newp->GetOutput());
	p->SetOutputParametersFromImage( rr->GetOutput() );
	p->Update();

	typename ImageWriter::Pointer w = ImageWriter::New();
	w->SetInput(p->GetOutput());

	ss.str("");
	ss << outPrefix << "_tpms";
	w->SetFileName(ss.str().c_str());
	w->Update();
}
