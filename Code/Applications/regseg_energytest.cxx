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

#include "regseg_energytest.hxx"

int main(int argc, char *argv[]) {
	std::string outfilename = "energies.json";
	std::string descfile, maskfile;
	std::vector< std::string > surfnames, refnames;
	std::stringstream ss;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("surfaces,S", bpo::value < std::vector< std::string > >(&surfnames)->multitoken()->required(), "input surfaces" )
			("reference,R", bpo::value< std::vector< std::string > >(&refnames)->multitoken()->required(), "reference image" )
			("descriptors,D", bpo::value< std::string >(&descfile), "descriptors file, if not present they will be computed from current segmentation" )
			("mask,M", bpo::value< std::string >(&maskfile), "inputmask" )
			("reorient", bpo::bool_switch(), "reorient images")
			("output,o", bpo::value < std::string > (&outfilename), "prefix for output files");

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

	// Some initial declarations
	bool reorient = vm.count("reorient") > 0;
	DirectionType idmat; idmat.SetIdentity();

	// Read surfaces
	InputMeshContainer priors;
	for( size_t i = 0; i<surfnames.size(); i++) {
		typename ReaderType::Pointer r = ReaderType::New();
		r->SetFileName(surfnames[i]);
		r->Update();
		priors.push_back(r->GetOutput());
	}

	// Read image channels and plug them into combine filter
	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();
	for (size_t i = 0; i < refnames.size(); i++ ) {
		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName( refnames[i] );
		r->Update();
		comb->SetInput(i,r->GetOutput());
	}
	comb->Update();
	ReferencePointer ref = comb->GetOutput();

	// Get some necessary features
	size_t ncomps = ref->GetNumberOfComponentsPerPixel();
	DirectionType native_dir = ref->GetDirection();
	PointType native_orig = ref->GetOrigin();

	// Reorient is necessary depending on the cli switch
	DirectionType reor_dir;
	PointType reor_orig;
	if (reorient) {
		typename Orienter::Pointer orient = Orienter::New();
		orient->UseImageDirectionOn();
		orient->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
		orient->SetInput(comb->GetOutput());
		orient->Update();
		ref = orient->GetOutput();

		reor_dir = ref->GetDirection();
		reor_orig = ref->GetOrigin();
		DirectionType itk; itk.SetIdentity();
		itk(0,0) = -1.0; itk(1,1) = -1.0;

		PointType neworig = itk * ref->GetOrigin();
		ref->SetDirection(idmat);
		ref->SetOrigin(neworig);
	}

	// Process mask (if available, generate empty otherwise)
	typename ChannelType::Pointer mask = ChannelType::New();
	mask->SetOrigin( ref->GetOrigin() );
	mask->SetDirection( ref->GetDirection() );
	mask->SetRegions( ref->GetLargestPossibleRegion().GetSize() );
	mask->SetSpacing( ref->GetSpacing() );
	mask->Allocate();
	mask->FillBuffer(1.0);

	if( vm.count("mask") ) {
		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName(maskfile);
		r->Update();
		typename ChannelType::Pointer inputmask = r->GetOutput();

		if (reorient) {
			typename ChannelOrienter::Pointer orient = ChannelOrienter::New();
			orient->UseImageDirectionOn();
			orient->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
			orient->SetInput(inputmask);
			orient->Update();
			inputmask = orient->GetOutput();
			inputmask->SetDirection(ref->GetDirection());
			inputmask->SetOrigin(ref->GetOrigin());
		}

		size_t npix = mask->GetLargestPossibleRegion().GetNumberOfPixels();
		ChannelPixelType* mbuff = mask->GetBufferPointer();
		const ChannelPixelType* msbuff = inputmask->GetBufferPointer();
		for (size_t pix = 0; pix < npix; pix++) {
			if ( *(msbuff+pix) > 1.0e-8 ) {
				*(mbuff+pix) = 0.0;
			}
		}
	}

	// Calculate sampling grid properties and initialize
	SizeType ref_size = ref->GetLargestPossibleRegion().GetSize();
	PointType ref_origin, ref_lastcenter, ref_end;
	ContinuousIndex tmp_idx;
	tmp_idx.Fill( -0.5 );
	ref->TransformContinuousIndexToPhysicalPoint( tmp_idx, ref_origin );
	for ( size_t dim = 0; dim<Dimension; dim++)  tmp_idx[dim]= ref_size[dim]-1.0;
	ref->TransformContinuousIndexToPhysicalPoint( tmp_idx, ref_lastcenter );
	for ( size_t dim = 0; dim<Dimension; dim++)  tmp_idx[dim]= ref_size[dim]- 0.5;
	ref->TransformContinuousIndexToPhysicalPoint( tmp_idx, ref_end );
	SizeType exp_size;
	for (size_t i = 0; i<Dimension; i++ ) exp_size[i] = (unsigned int) (ref_size[i] * 2);

	PointType firstcenter;
	SpacingType spacing;
	VectorType step;

	for (size_t i = 0; i<Dimension; i++ ){
		step[i] = (ref_end[i] - ref_origin[i]) / (1.0*exp_size[i]);
		spacing[i]= fabs( step[i] );
		firstcenter[i] = ref_origin[i] + 0.5 * step[i];
	}

	ReferencePointer refSamplingGrid = ReferenceImageType::New();
	refSamplingGrid->SetNumberOfComponentsPerPixel(ncomps);
	refSamplingGrid->SetOrigin( firstcenter );
	refSamplingGrid->SetDirection( idmat );
	refSamplingGrid->SetRegions( exp_size );
	refSamplingGrid->SetSpacing( spacing );
	refSamplingGrid->Allocate();

	// Binarization & downsamplings
	BinarizeMeshFilterPointer newp = BinarizeMeshFilterType::New();
	newp->SetInputs( priors );
	newp->SetOutputReference( refSamplingGrid );
	newp->Update();

	DownsamplePointer p = DownsampleFilter::New();
	p->SetInput(newp->GetOutput());
	p->SetOutputParametersFromImage( ref );
	p->Update();

	typename ProbmapType::Pointer tpms = p->GetOutput();

	if( reorient ) {
		tpms->SetDirection(reor_dir);
		tpms->SetOrigin(reor_orig);
		typename ProbmapsOrienter::Pointer orient_tpm = ProbmapsOrienter::New();
		orient_tpm->UseImageDirectionOn();
		orient_tpm->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
		orient_tpm->SetInput(tpms);
		orient_tpm->Update();

		tpms = orient_tpm->GetOutput();
		tpms->SetDirection(native_dir);
		tpms->SetOrigin(native_orig);
	}


	typename ImageWriter::Pointer w = ImageWriter::New();
	w->SetInput(tpms);
	w->SetFileName("tpms");
	w->Update();

	EnergyModelPointer model = EnergyModelType::New();
	if (vm.count("descriptors")) {
		model->ReadDescriptorsFromFile(descfile);
	} else {
		model->SetInput(ref);
		model->SetPriorsMap(p->GetOutput());
		model->SetMask(mask);
		model->Update();

		std::string jsonstr = model->PrintFormattedDescriptors();
		std::ofstream outdescfile("descriptors.json");
		outdescfile << jsonstr;
		outdescfile.close();
	}

	EnergyFilterPointer ecalc = EnergyFilter::New();
	ecalc->SetInput(ref);
	ecalc->SetPriorsMap(p->GetOutput());
	ecalc->SetMask(mask);
	ecalc->SetModel(model);
	ecalc->Update();
	MeasureArray values = ecalc->GetEnergies();
	MeasureType total = 0.0;

	Json::Value root = Json::Value( Json::objectValue );
	root["regions"] = Json::Value( Json::arrayValue );

	for( size_t roi = 0; roi < values.Size(); roi++ ) {
		root["regions"].append(Json::Value(values[roi]));
		total+= values[roi];
	}

	root["total"]= Json::Value(total);
	std::ofstream outfile(outfilename.c_str());
	outfile << root.toStyledString();
	outfile.close();
	return EXIT_SUCCESS;
}
