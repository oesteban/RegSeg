/*
 Copyright (c) 2012, Oscar Esteban - oesteban@die.upm.es
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
	std::string logFileName = "regseg.log";
	bool outImages = false;
	size_t verbosity = 1;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
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


	// Set up log file
	std::ofstream logfile((outPrefix + logFileName ).c_str());

	// Read fixed image(s) --------------------------------------------------------------
	clock_t preProcessStart = clock();


	// Initialize LevelSet function
	FunctionalType::Pointer functional = FunctionalType::New();
	// Connect Optimizer
	OptimizerPointer opt = Optimizer::New();
	opt->SetFunctional( functional );


	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();
	logfile << " * Target feature, number of components= " << fixedImageNames.size() << "." << std::endl;
	for (size_t i = 0; i < fixedImageNames.size(); i++ ) {
		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();
		comb->SetInput(i,r->GetOutput());

		logfile << "\t* Input " << i << ": " << fixedImageNames[i] << "." << std::endl;
	}

	ChannelType::DirectionType dir; dir.SetIdentity();
	comb->Update();
	ImageType::Pointer im = comb->GetOutput();
	im->SetDirection( dir );
	functional->SetReferenceImage( im );


	// Read moving surface(s) -----------------------------------------------------------
	logfile << " * Number of moving surfaces= " << movingSurfaceNames.size() << "." << std::endl;
	for (size_t i = 0; i < movingSurfaceNames.size(); i++) {
		ReaderType::Pointer polyDataReader = ReaderType::New();
		polyDataReader->SetFileName( movingSurfaceNames[i] );
		polyDataReader->Update();
		functional->AddShapePrior( polyDataReader->GetOutput() );

		logfile << "\t* Mesh " << i << ": " << movingSurfaceNames[i] << "." << std::endl;
	}

	// Set up registration ------------------------------------------------------------
	if (vmap.count("grid-size,g")) {
		opt->SetGridSize( vmap["grid-size,g"].as<int>() );
	}
	if (vmap.count("iterations,i")) {
		opt->SetNumberOfIterations( vmap["iterations,i"].as<int>() );
	}
	if (vmap.count("step-size,s")) {
		opt->SetStepSize( vmap["step-size,s"].as<float>() );

	}
	if (vmap.count("alpha,a")) {
		opt->SetAlpha( vmap["alpha,a"].as<float>() );
	}

	clock_t preProcessStop = clock();
	float pre_tot_t = (float) (((double) (preProcessStop - preProcessStart))
			/ CLOCKS_PER_SEC);
	int h, min, sec;
	h = (pre_tot_t / 3600);
	min = (((int) pre_tot_t) % 3600) / 60;
	sec = (((int) pre_tot_t) % 3600) % 60;

	char pre_time[50];
	sprintf(pre_time, "\t* Pre-processing Total Time = %02d:%02d:%02d hours\n",
			h, min, sec);


	// Start registration -------------------------------------------------------------
	logfile << " --------------------------------- Starting registration process." << std::endl;
	clock_t initTime = clock();

	opt->Start();

	clock_t finishTime = clock();
	logfile << " --------------------------------- Finished registration process." << std::endl;

	logfile << "\n Summary:" << std::endl;

	logfile << "\t* Final energy = " << functional->GetValue() << "." << std::endl;


	float tot_t = (float) (((double) (finishTime - initTime)) / CLOCKS_PER_SEC);
	h = (tot_t / 3600);
	min = (((int) tot_t) % 3600) / 60;
	sec = (((int) tot_t) % 3600) % 60;

	char reg_time[50];
	sprintf(reg_time, "\t* Registration Total Time = %02d:%02d:%02d hours\n", h,
			min, sec);

	logfile << pre_time << std::endl;
	logfile << reg_time << std::endl;


	// Write out final results ---------------------------------------------------------
    size_t nCont = functional->GetCurrentContourPosition().size();
    for ( size_t contid = 0; contid < nCont; contid++) {
        std::stringstream ss;
        ss << "final-cont0" << contid << ".vtk";
    	WriterType::Pointer polyDataWriter = WriterType::New();
    	polyDataWriter->SetInput( functional->GetCurrentContourPosition()[contid] );
    	polyDataWriter->SetFileName( "deformed2-wm.vtk" );
    	polyDataWriter->Update();
    }


	return EXIT_SUCCESS;
}
