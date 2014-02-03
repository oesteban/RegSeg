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

#include <boost/shared_ptr.hpp>

int main(int argc, char *argv[]) {
	std::string outPrefix;
	std::vector< std::string > fixedImageNames, movingSurfaceNames;
	std::string logFileName = ".log";
	bool outImages = false;
	size_t verbosity = 1;

	bpo::options_description all_desc("Usage");
	bpo::options_description general_desc("General options");
	general_desc.add_options()
			("help,h", "show help message")
			("fixed-images,F", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("moving-surfaces,M", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken()->required(),	"moving image file")
			("transform-levels,L", bpo::value< size_t > (), "number of multi-resolution levels for the transform")
			("output-prefix,o", bpo::value < std::string > (&outPrefix), "prefix for output files")
			("output-all", bpo::bool_switch(&outImages),"output intermediate images")
			("logfile,l", bpo::value<std::string>(&logFileName), "log filename")
			("verbosity,V", bpo::value<size_t>(&verbosity), "verbosity level ( 0 = no output; 5 = verbose )");

	bpo::options_description level_desc("Registration level options");
	add_level_options( level_desc );


	std::vector< std::string > cli_token;
	std::vector< std::string > cli_general;
	std::vector< std::vector< std::string > > cli_levels;

	bool isLevel = false;
	std::string token;
	for( int i = 0; i < argc; i++) {
		token = argv[i];
		if( !isLevel ) {
			if ( token.at(0)=='[' ) {
				cli_token.clear();
				isLevel = true;
				token.erase(0,1);
			} else {
				cli_general.push_back( argv[i] );
			}
		}

		if( isLevel ) {
			if ( *token.rbegin() == ']' ) {
				token = token.substr( 0, token.size() -1 );
				isLevel = false;
			}
			cli_token.push_back( token );

			if ( !isLevel ) {
				cli_levels.push_back( cli_token );
			}
		}
	}

	size_t cli_nlevels = cli_levels.size();

	bpo::variables_map vm_general;
	std::vector< bpo::variables_map > vm_levels;


	all_desc.add( general_desc ).add( level_desc );

	try {
		// Deal with general options
		if( cli_nlevels == 0 ) {
			bpo::store(	bpo::command_line_parser( cli_general ).options(all_desc).run(),vm_general);
		} else {
			bpo::store( bpo::command_line_parser( cli_general ).options(general_desc).run(), vm_general );
		}


		if (vm_general.count("help")) {
			std::cout << all_desc << std::endl;
			return 1;
		}

		bpo::notify(vm_general);


		if( cli_nlevels > 0 ) {
			for ( size_t i = 0; i<cli_nlevels; i++ ) {
				bpo::variables_map vm;
				bpo::options_description ndesc("Level " + boost::lexical_cast<std::string> (i) + " options");
				add_level_options( ndesc );

				bpo::store(	bpo::command_line_parser( cli_levels[i] ).options(ndesc).run(), vm );
				bpo::notify( vm );
				vm_levels.push_back( vm );
			}
		}
	} catch (boost::exception_detail::clone_impl
			< boost::exception_detail::error_info_injector<	boost::program_options::required_option> > &err) {
		std::cout << "Error: " << err.what() << std::endl;
		std::cout << all_desc << std::endl;
		return EXIT_FAILURE;
	}

	// Create the JSON output object
	Json::Value root;
	root["description"] = "RSTK Summary File";
	std::time_t time;
	root["information"] = std::ctime( &time );

	// Initialize registration
	RegistrationPointer acwereg = RegistrationType::New();


	// Read target feature(s) -----------------------------------------------------------
	root["inputs"]["target"]["components"]["size"] = Json::Int (fixedImageNames.size());
	root["inputs"]["target"]["components"]["type"] = std::string("feature");
	Json::Value targetjson(Json::arrayValue);

	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();
	for (size_t i = 0; i < fixedImageNames.size(); i++ ) {
		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();
		comb->SetInput(i,r->GetOutput());
		targetjson.append( fixedImageNames[i] );
	}
	root["inputs"]["target"]["components"] = targetjson;

    comb->Update();
	acwereg->SetFixedImage( comb->GetOutput() );


	// Read moving surface(s) -----------------------------------------------------------
	root["inputs"]["moving"]["components"]["size"] = Json::Int (movingSurfaceNames.size());
	root["inputs"]["moving"]["components"]["type"] = std::string("surface");
	Json::Value movingjson(Json::arrayValue);

	for (size_t i = 0; i < movingSurfaceNames.size(); i++) {
		ReaderType::Pointer polyDataReader = ReaderType::New();
		polyDataReader->SetFileName( movingSurfaceNames[i] );
		polyDataReader->Update();
		acwereg->AddShapePrior( polyDataReader->GetOutput() );
		movingjson.append( movingSurfaceNames[i] );
	}
	root["inputs"]["moving"]["components"] = movingjson;

	// Set up registration ------------------------------------------------------------
	if ( vm_general.count("transform-levels") && cli_nlevels == 0 ) {
		acwereg->SetNumberOfLevels( vm_general["transform-levels"].as<size_t>() );
		acwereg->SetUseGridLevelsInitialization( true );
	}
	if ( cli_nlevels ) {
		acwereg->SetNumberOfLevels( cli_nlevels );
		acwereg->SetUseGridLevelsInitialization( true );
	}

	for( size_t i = 0; i < cli_nlevels; i++ ) {
		bpo::variables_map vm = vm_levels[i];

		if (vm.count("iterations")) {
			size_t nit = vm["iterations"].as< size_t >();
			acwereg->SetNumberOfIterationsElement( i, nit );
		}
        //
		//if (vm.count("step-size")) {
		//	std::vector< double > ssize = vm["step-size"].as< std::vector< double > >();
		//	acwereg->SetStepSizeElement( i,  ssize[0] );
        //
		//}
		//if (vm.count("alpha")) {
		//	acwereg->SetAlphaElement( i, vm["alpha"].as<float>() );
		//}
        //
		//if (vm.count("beta")) {
		//	acwereg->SetBetaElement( i, vm["beta"].as<float>() );
		//}
		//if (vm.count("descriptors-update-iterations")) {
		//	size_t updDesc =  vm["descriptors-update-iterations"].as<size_t>();
		//	acwereg->FillDescriptorRecomputationFreq(updDesc);
		//}
		//if (vm.count("grid-size")) {
		//	acwereg->SetGridSize( vm["grid-size"].as<size_t>() );
		//}
	}

	acwereg->Update();

	//
	// Write out final results ---------------------------------------------------------
	//

	// Displacementfield
	//typename DisplacementFieldWriter::Pointer p = DisplacementFieldWriter::New();
	//p->SetFileName( (outPrefix + "_field.nii.gz" ).c_str() );
	//p->SetInput( opt->GetCurrentDisplacementField() );
	//p->Update();
    //
	//// Contours and regions
    //size_t nCont = functional->GetCurrentContours().size();
    //for ( size_t contid = 0; contid < nCont; contid++) {
    //	bfs::path contPath(movingSurfaceNames[contid]);
    //	WriterType::Pointer polyDataWriter = WriterType::New();
    //	polyDataWriter->SetInput( functional->GetCurrentContours()[contid] );
    //	polyDataWriter->SetFileName( (outPrefix + "_" + contPath.filename().string()).c_str() );
    //	polyDataWriter->Update();
    //
    //	typename ROIWriter::Pointer w = ROIWriter::New();
    //	w->SetInput( functional->GetCurrentRegion(contid) );
    //	w->SetFileName( (outPrefix + "_roi_" + contPath.stem().string() + ".nii.gz" ).c_str() );
    //	w->Update();
    //}
    //
    //// Last ROI (excluded region)
	//typename ROIWriter::Pointer w = ROIWriter::New();
	//w->SetInput( functional->GetCurrentRegion(nCont) );
	//w->SetFileName( (outPrefix + "_roi_background.nii.gz" ).c_str() );
	//w->Update();


	// Set-up & write out log file
	std::ofstream logfile((outPrefix + logFileName ).c_str());
	logfile << root;

	return EXIT_SUCCESS;
}
