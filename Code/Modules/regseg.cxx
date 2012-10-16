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
	bool outImages = false;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
			("fixed-images,F", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("moving-surfaces,M", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken()->required(),	"moving image file")
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("output-all", bpo::bool_switch(&outImages),"output intermediate images");
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
	std::ofstream logfile((outPrefix + "regseg.log").c_str());

	// Read fixed image(s) --------------------------------------------------------------
	clock_t preProcessStart = clock();

	ImageVectorType::ConstPointer fixedImages;
	if ( fixedImageNames.size() == 1 ) {
		ImageReader::Pointer fixedImageReader = ImageReader::New();
		fixedImageReader->SetFileName(fixedImageNames[0]);
		fixedImageReader->Update();
		fixedImages = fixedImageReader->GetOutput();
		logfile << " * Loaded unique input " << ": " << fixedImageNames[0] << "." << std::endl;
	} else {
		std::vector< const ImageComponentType* > fixedImagesVector;

		ImageToVectorFilterType::Pointer inputVectorFilter = ImageToVectorFilterType::New();
		for (size_t i = 0; i < fixedImageNames.size(); i++ ) {
			ImageComponentReader::Pointer fixedImageReader = ImageComponentReader::New();
			fixedImageReader->SetFileName(fixedImageNames[i]);
			fixedImageReader->Update();
			//fixedImagesVector.push_back( fixedImageReader->GetOutput() );
			inputVectorFilter->SetNthInput(i, fixedImageReader->GetOutput() );
			inputVectorFilter->Update();
			logfile << " * Loaded input " << i << ": " << fixedImageNames[i] << "." << std::endl;
		}
		fixedImages = inputVectorFilter->GetOutput();

		logfile << " * Number of Components= " << fixedImages->GetNumberOfComponentsPerPixel() << "." << std::endl;
	}

	// QUESTION Read fixed image parameters -------------------------------------------------


	// Read moving surface(s) -----------------------------------------------------------
	SurfaceVector movingSurfaces;
	for (size_t i = 0; i < movingSurfaceNames.size(); i++) {
		SurfaceReader::Pointer surfReader = SurfaceReader::New();
		surfReader->SetFileName( movingSurfaceNames[i] );
		surfReader->Update();
		movingSurfaces.push_back( surfReader->GetOutput() );
	}


	// Set up registration ------------------------------------------------------------

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
	clock_t initTime = clock();
	// ...

	clock_t finishTime = clock();
	float tot_t = (float) (((double) (finishTime - initTime)) / CLOCKS_PER_SEC);
	h = (tot_t / 3600);
	min = (((int) tot_t) % 3600) / 60;
	sec = (((int) tot_t) % 3600) % 60;

	char reg_time[50];
	sprintf(reg_time, "\t* Registration Total Time = %02d:%02d:%02d hours\n", h,
			min, sec);

	logfile << pre_time << std::endl;
	logfile << reg_time << std::endl;

	return EXIT_SUCCESS;
}
