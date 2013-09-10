/*
 Copyright (c) 2012, Oscar Esteban - oesteban@dionte.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM)
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the names of the BIT-UPM and LTS5-EPFL, nor the names of its contributors
 may be used to endorse or promote products derived from this software
 without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "regseg.h"

int main(int argc, char *argv[]) {
	std::string outPrefix;
	std::vector< std::string > fixedImageNames, movingSurfaceNames;
	std::string logFileName = ".log";
	bool outImages = false;
	size_t verbosity = 1;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help,h", "show help message")
			("fixed-images,F", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("moving-surfaces,M", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken()->required(),	"moving image file")
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("output-all", bpo::bool_switch(&outImages),"output intermediate images")
			("alpha,a", bpo::value< float > (), "alpha value in regularization")
			("beta,b", bpo::value< float > (), "beta value in regularization")
			("step-size,s", bpo::value< float > (), "step-size value in optimization")
			("iterations,i", bpo::value< int > (), "number of iterations")
			("grid-size,g", bpo::value< int > (), "grid size")
			("logfile,l", bpo::value<std::string>(&logFileName), "log filename")
			("verbosity,V", bpo::value<size_t>(&verbosity), "verbosity level ( 0 = no output; 5 = verbose )");
	bpo::positional_options_description pdesc;
	bpo::variables_map vmap;
	bpo::store(
			bpo::command_line_parser(argc, argv).options(desc).positional(pdesc).run(),
			vmap);

	if (vmap.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}

	try {
		bpo::notify(vmap);
	} catch (boost::exception_detail::clone_impl
			< boost::exception_detail::error_info_injector<	boost::program_options::required_option> > &err) {
		std::cout << "Error: " << err.what() << std::endl;
		std::cout << desc << std::endl;
		return EXIT_FAILURE;
	}

	// Create the JSON output object
	Json::Value root;
	root["description"] = "RSTK Summary File";
	std::time_t time;
	root["information"] = std::ctime( &time );

	// Read fixed image(s) --------------------------------------------------------------
	clock_t preProcessStart = clock();


	// Initialize LevelSet function
	FunctionalType::Pointer functional = FunctionalType::New();
	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetFunctional( functional );


	//typename ObserverType::Pointer o = ObserverType::New();
	//o->SetCallbackFunction( opt, & PrintIteration );
	//o->AddObserver( itk::IterationEvent(), o );

	// Read target feature(s) -----------------------------------------------------------
	root["inputs"]["target"]["components"]["size"] = Json::Int (fixedImageNames.size());
	root["inputs"]["target"]["components"]["type"] = std::string("feature");
	Json::Value targetjson(Json::arrayValue);

	ImageType::DirectionType dir; dir.SetIdentity();
	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();
	for (size_t i = 0; i < fixedImageNames.size(); i++ ) {
		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();
		comb->SetInput(i,r->GetOutput());

		targetjson.append( fixedImageNames[i]);
	}
	root["inputs"]["target"]["components"] = targetjson;

	ImageType::Pointer im = comb->GetOutput();
	functional->SetReferenceImage( im );

	// Read moving surface(s) -----------------------------------------------------------
	root["inputs"]["moving"]["components"]["size"] = Json::Int (movingSurfaceNames.size());
	root["inputs"]["moving"]["components"]["type"] = std::string("surface");
	Json::Value movingjson(Json::arrayValue);

	for (size_t i = 0; i < movingSurfaceNames.size(); i++) {
		ReaderType::Pointer polyDataReader = ReaderType::New();
		polyDataReader->SetFileName( movingSurfaceNames[i] );
		polyDataReader->Update();
		functional->AddShapePrior( polyDataReader->GetOutput() );
		movingjson.append( movingSurfaceNames[i] );
	}
	root["inputs"]["moving"]["components"] = movingjson;

	// Set up registration ------------------------------------------------------------
	if (vmap.count("grid-size")) {
		opt->SetGridSize( vmap["grid-size"].as<int>() );
	}
	if (vmap.count("iterations")) {
		opt->SetNumberOfIterations( vmap["iterations"].as<int>() );
	}
	if (vmap.count("step-size")) {
		opt->SetStepSize( vmap["step-size"].as<float>() );

	}
	if (vmap.count("alpha")) {
		opt->SetAlpha( vmap["alpha"].as<float>() );
	}

	clock_t preProcessStop = clock();
	float pre_tot_t = (float) (((double) (preProcessStop - preProcessStart)) / CLOCKS_PER_SEC);
	root["time"]["preprocessing"] = pre_tot_t;

	// Start registration -------------------------------------------------------------
	std::cout << " --------------------------------- Starting registration process." << std::endl;
	clock_t initTime = clock();

	opt->Start();

	clock_t finishTime = clock();
	std::cout << " --------------------------------- Finished registration process." << std::endl;

	float tot_t = (float) (((double) (finishTime - initTime)) / CLOCKS_PER_SEC);
	root["time"]["processing"] = tot_t;

	//
	// Write out final results ---------------------------------------------------------
	//

	// Displacementfield
	typename DisplacementFieldWriter::Pointer p = DisplacementFieldWriter::New();
	p->SetFileName( (outPrefix + "_field.nii.gz" ).c_str() );
	p->SetInput( functional->GetCurrentDisplacementField() );
	p->Update();

	// Contours and regions
    size_t nCont = functional->GetCurrentContourPosition().size();
    for ( size_t contid = 0; contid < nCont; contid++) {
    	bfs::path contPath(movingSurfaceNames[contid]);
    	WriterType::Pointer polyDataWriter = WriterType::New();
    	polyDataWriter->SetInput( functional->GetCurrentContourPosition()[contid] );
    	polyDataWriter->SetFileName( (outPrefix + "_" + contPath.filename().string()).c_str() );
    	polyDataWriter->Update();

    	typename ROIWriter::Pointer w = ROIWriter::New();
    	w->SetInput( functional->GetCurrentRegion(contid) );
    	w->SetFileName( (outPrefix + "_roi_" + contPath.stem().string() + ".nii.gz" ).c_str() );
    	w->Update();
    }

    // Last ROI (excluded region)
	typename ROIWriter::Pointer w = ROIWriter::New();
	w->SetInput( functional->GetCurrentRegion(nCont) );
	w->SetFileName( (outPrefix + "_roi_background.nii.gz" ).c_str() );
	w->Update();


	// JSON Summary
	root["summary"]["energy"]["total"] = opt->GetCurrentMetricValue();
	root["summary"]["energy"]["data"] = functional->GetValue();
	root["summary"]["energy"]["regularization"] = opt->GetCurrentRegularizationEnergy();
	root["summary"]["iterations"] = Json::Int (opt->GetCurrentIteration());
	root["summary"]["conv_status"] = opt->GetStopCondition();
	root["summary"]["stop_msg"] = opt->GetStopConditionDescription();
	Json::Value dumbArray(Json::arrayValue);
	root["summary"]["energy"]["evolution_tot"] = dumbArray;
	root["summary"]["energy"]["evolution_dat"] = dumbArray;
	root["summary"]["energy"]["evolution_reg"] = dumbArray;

	// Set-up & write out log file
	std::ofstream logfile((outPrefix + logFileName ).c_str());
	logfile << root;

	return EXIT_SUCCESS;
}
