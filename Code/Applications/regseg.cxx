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
#include <boost/asio/signal_set.hpp>


int main(int argc, char *argv[]) {
	std::string outPrefix;
	std::vector< std::string > fixedImageNames, movingSurfaceNames, targetSurfaceNames;
	std::string logFileName = "";
	bool outImages = false;

	bpo::options_description all_desc("Usage");
	bpo::options_description general_desc("General options");
	general_desc.add_options()
			("help,h", "show help message")
			("fixed-images,F", bpo::value < std::vector<std::string>	> (&fixedImageNames)->multitoken()->required(), "fixed image file")
			("surface-priors,P", bpo::value < std::vector<std::string>	> (&movingSurfaceNames)->multitoken()->required(),	"shape priors")
			("surface-target,T", bpo::value < std::vector<std::string>	> (&targetSurfaceNames)->multitoken(),	"final shapes to evaluate metric (only testing purposes)")
			("fixed-mask,M", bpo::value< std::string >(), "fixed image mask")
			("transform-levels,L", bpo::value< size_t > (), "number of multi-resolution levels for the transform")
			("output-prefix,o", bpo::value < std::string > (&outPrefix)->default_value("regseg"), "prefix for output files")
			("logfile,l", bpo::value<std::string>(&logFileName), "log filename")
			("monitoring-verbosity,v", bpo::value<size_t>()->default_value(DEFAULT_VERBOSITY), "verbosity level of intermediate results monitoring ( 0 = no output; 5 = verbose )");

	bpo::options_description opt_desc("Optimizer options (by levels)");
	OptimizerType::AddOptions( opt_desc );
	bpo::options_description fun_desc("Functional options (by levels)");
	FunctionalType::AddOptions( fun_desc );


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


	all_desc.add( general_desc ).add( opt_desc ).add( fun_desc );

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
				OptimizerType::AddOptions( ndesc );
				FunctionalType::AddOptions( ndesc );

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

	// Initialize registration
	RegistrationPointer acwereg = RegistrationType::New();
	acwereg->SetOutputPrefix( outPrefix );
	acwereg->SetVerbosity( vm_general["monitoring-verbosity"].as< size_t >() );

	// Create the JSON output object
	Json::Value root;
	root["description"]["title"] = "ACWE-Reg Summary File";
	std::time_t rawtime;
	std::tm* timeinfo;
	char buffer[40];
	std::time(&rawtime);
	timeinfo = std::localtime(&rawtime);
	std::strftime( buffer, 20, "%Y-%m-%d",timeinfo);
	root["description"]["date"] = buffer;
	std::strftime( buffer, 20, "%H:%M:%S",timeinfo);
	root["description"]["time"] = buffer;



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

	// Read target mask -----------------------------------------------------------------
	if( vm_general.count( "fixed-mask" ) ) {
		std::string maskfname = vm_general["fixed-mask"].as< std::string >();
		root["inputs"]["target"]["mask"] = maskfname;

		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName(maskfname);
		r->Update();
		acwereg->SetFixedMask( r->GetOutput() );
	}


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

	for (size_t i = 0; i < targetSurfaceNames.size(); i++) {
		ReaderType::Pointer polyDataReader = ReaderType::New();
		polyDataReader->SetFileName( targetSurfaceNames[i] );
		polyDataReader->Update();
		acwereg->AddShapeTarget( polyDataReader->GetOutput() );
	}

	// Set up registration ------------------------------------------------------------
	if ( vm_general.count("transform-levels") && cli_nlevels == 0 ) {
		acwereg->SetNumberOfLevels( vm_general["transform-levels"].as<size_t>() );
		acwereg->SetUseGridLevelsInitialization( true );
	}
	if ( cli_nlevels ) {
		acwereg->SetNumberOfLevels( cli_nlevels );
		acwereg->SetUseGridLevelsInitialization( false );
	}

	for( size_t i = 0; i < cli_nlevels; i++ ) {
		bpo::variables_map vm = vm_levels[i];
		acwereg->SetSettingsOfLevel( i, vm );
	}

	try {
		acwereg->Update();
	} catch (const std::exception &exc) {
		// Set-up & write out log file
		root["levels"] = acwereg->GetJSONRoot();
		root["error"]["message"] = std::string(exc.what());
		std::ofstream logfile((outPrefix + logFileName + ".log" ).c_str());
		logfile << root;
		throw;
	} catch (...) {
		root["levels"] = acwereg->GetJSONRoot();
		root["error"]["message"] = std::string("Unknown exception");
		std::ofstream logfile((outPrefix + logFileName + ".log" ).c_str());
		logfile << root;
		throw;
	}

	root["levels"] = acwereg->GetJSONRoot();
	// Set-up & write out log file
	std::ofstream logfile((outPrefix + logFileName + ".log" ).c_str());
	logfile << root;

	//
	// Write out final results ---------------------------------------------------------
	//
	size_t nlevels = acwereg->GetNumberOfLevels();

	for (size_t i = 0; i < nlevels; i++) {
		std::stringstream sb;
		sb << outPrefix << "_coeff_" << i << ".vtu";
		typename CoeffWriter::Pointer w = CoeffWriter::New();
		w->SetFileName(sb.str().c_str());
		w->SetInput(acwereg->GetOptimizerOfLevel(i)->GetTransform()->GetFlatParameters());
		w->Update();
	}

	// Displacementfield
	typename FieldWriter::Pointer fwrite = FieldWriter::New();
	fwrite->SetFileName( (outPrefix + "_field.nii.gz" ).c_str() );
	fwrite->SetInput( acwereg->GetDisplacementField() );
	fwrite->Update();

	// Contours and regions
	ContourList conts = acwereg->GetCurrentContours();
    size_t nCont = conts.size();
    for ( size_t contid = 0; contid < nCont; contid++) {
    	bfs::path contPath(movingSurfaceNames[contid]);
    	WriterType::Pointer polyDataWriter = WriterType::New();
    	std::stringstream ss;
    	ss << outPrefix << "_swarped_" << contid << ".vtk";
    	polyDataWriter->SetInput( conts[contid] );
    	polyDataWriter->SetFileName( ss.str().c_str() );
    	polyDataWriter->Update();
    }

    //for ( size_t contid = 0; contid <= nCont; contid++) {
	//	typename ProbabilityMapWriter::Pointer w = ProbabilityMapWriter::New();
	//	w->SetInput( acwereg->GetFunctionalOfLevel(-1)->GetCurrentMap(contid) );
	//	std::stringstream ss;
	//	ss << outPrefix << "_final_tpm_" << contid << ".nii.gz";
	//	w->SetFileName( ss.str().c_str() );
	//	w->Update();
    //}

	// Read and transform images if present
	for( size_t i = 0; i<fixedImageNames.size(); i++) {
		typename ImageReader::Pointer r = ImageReader::New();
		r->SetFileName( fixedImageNames[i] );
		r->Update();

		typename ChannelType::Pointer im = r->GetOutput();
		typename ChannelType::DirectionType dir = im->GetDirection();
		typename ChannelType::PointType ref_orig = im->GetOrigin();

		typename ChannelType::DirectionType itk;
		itk.SetIdentity();
		itk(0,0)=-1.0;
		itk(1,1)=-1.0;
		im->SetDirection( dir * itk );
		im->SetOrigin( itk * ref_orig );

		WarpFilterPointer res = WarpFilter::New();
		res->SetInterpolator( DefaultInterpolator::New() );
		res->SetOutputParametersFromImage( im );
		res->SetInput( im );
		res->SetDisplacementField( acwereg->GetDisplacementField() );
		res->Update();


		ThresholdPointer th = ThresholdFilter::New();
		th->SetInput( res->GetOutput() );
		th->ThresholdBelow( 0.0 );
		th->SetOutsideValue( 0.0 );

		typename ChannelType::Pointer im_res = th->GetOutput();
		im_res->SetDirection( dir );
		im_res->SetOrigin( ref_orig );

		std::stringstream ss;
		ss.str("");
		ss << outPrefix << "_warped_" << i << ".nii.gz";
		typename ImageWriter::Pointer w = ImageWriter::New();
		w->SetInput( im_res );
		w->SetFileName( ss.str().c_str() );
		w->Update();

	}

	return EXIT_SUCCESS;
}
