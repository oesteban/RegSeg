// --------------------------------------------------------------------------------------
// File:          regseg_energytest.cxx
// Date:          Feb 12, 2015
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "regseg_energytest.hxx"

int main(int argc, char *argv[]) {
	std::string outPrefix = "surf2vol";
	std::string descfile;
	std::vector< std::string > surfnames, refnames;

	bpo::options_description all_desc("Usage");
	all_desc.add_options()
			("help,h", "show help message")
			("surfaces,S", bpo::value < std::vector< std::string > >(&surfnames)->multitoken()->required(), "input surfaces" )
			("reference,R", bpo::value< std::vector< std::string > >(&refnames)->required(), "reference image" )
			("descriptors,D", bpo::value< std::string >(&descfile), "descriptors file, if not present they will be computed from current segmentation" )
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

	EnergyModelPointer model = EnergyModelType::New();
	if (vm.count("descriptors")) {
		model->ReadDescriptorsFromFile(descfile);
	}

	InputMeshContainer priors;

	for( size_t i = 0; i<surfnames.size(); i++) {
		typename ReaderType::Pointer r = ReaderType::New();
		r->SetFileName(surfnames[i]);
		r->Update();
		priors.push_back(r->GetOutput());
	}

	InputToVectorFilterType::Pointer comb = InputToVectorFilterType::New();
	for (size_t i = 0; i < refnames.size(); i++ ) {
		ImageReader::Pointer r = ImageReader::New();
		r->SetFileName( refnames[i] );
		r->Update();
		comb->SetInput(i,r->GetOutput());
	}
	comb->Update();
	size_t ncomps = comb->GetOutput()->GetNumberOfComponentsPerPixel();

	typename Orienter::Pointer orient = Orienter::New();
	orient->UseImageDirectionOn();
	orient->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
	orient->SetInput(comb->GetOutput());
	orient->Update();
	ReferencePointer ref = orient->GetOutput();
	SizeType ref_size = ref->GetLargestPossibleRegion().GetSize();

	DirectionType reor_dir = ref->GetDirection();
	PointType reor_orig = ref->GetOrigin();


	DirectionType idmat; idmat.SetIdentity();
	DirectionType itk; itk.SetIdentity();
	itk(0,0) = -1.0; itk(1,1) = -1.0;

	PointType neworig = itk * ref->GetOrigin();
	ref->SetDirection(idmat);
	ref->SetOrigin(neworig);

	PointType ref_origin, ref_lastcenter, ref_end;

	ContinuousIndex tmp_idx;
	tmp_idx.Fill( -0.5 );
	ref->TransformContinuousIndexToPhysicalPoint( tmp_idx, ref_origin );

	for ( size_t dim = 0; dim<Dimension; dim++)  tmp_idx[dim]= ref_size[dim]-1.0;
	ref->TransformContinuousIndexToPhysicalPoint( tmp_idx, ref_lastcenter );

	for ( size_t dim = 0; dim<Dimension; dim++)  tmp_idx[dim]= ref_size[dim]- 0.5;
	ref->TransformContinuousIndexToPhysicalPoint( tmp_idx, ref_end );

	SizeType exp_size;
	for (size_t i = 0; i<Dimension; i++ ){
		exp_size[i] = (unsigned int) (ref_size[i] * 2);
	}

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

	BinarizeMeshFilterPointer newp = BinarizeMeshFilterType::New();
	newp->SetInputs( priors );
	newp->SetOutputReference( refSamplingGrid );
	newp->Update();

	typename SegmentationType::Pointer seg = newp->GetOutputSegmentation();
	seg->SetDirection(reor_dir);
	seg->SetOrigin(reor_orig);

	typename SegmentationOrienter::Pointer orient_seg = SegmentationOrienter::New();
	orient_seg->UseImageDirectionOn();
	orient_seg->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
	orient_seg->SetInput(seg);
	orient_seg->Update();
	typename SegmentationType::Pointer seg_oriented = orient_seg->GetOutput();
	seg_oriented->SetDirection(comb->GetOutput()->GetDirection());
	seg_oriented->SetOrigin(comb->GetOutput()->GetOrigin());

	std::stringstream ss;
	ss << outPrefix << "_seg.nii.gz";
	typename SegmentationWriter::Pointer sw = SegmentationWriter::New();
	sw->SetFileName(ss.str().c_str());
	sw->SetInput(seg_oriented);
	sw->Update();

	DownsamplePointer p = DownsampleFilter::New();
	p->SetInput(newp->GetOutput());
	p->SetOutputParametersFromImage( ref );
	p->Update();

	typename ProbmapType::Pointer rawtpms = p->GetOutput();
	rawtpms->SetDirection(reor_dir);
	rawtpms->SetOrigin(reor_orig);

	typename ProbmapsOrienter::Pointer orient_tpm = ProbmapsOrienter::New();
	orient_tpm->UseImageDirectionOn();
	orient_tpm->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
	orient_tpm->SetInput(rawtpms);
	orient_tpm->Update();
	typename ProbmapType::Pointer tpm_oriented = orient_tpm->GetOutput();
	tpm_oriented->SetDirection(comb->GetOutput()->GetDirection());
	tpm_oriented->SetOrigin(comb->GetOutput()->GetOrigin());

	typename ImageWriter::Pointer w = ImageWriter::New();
	w->SetInput(tpm_oriented);

	ss.str("");
	ss << outPrefix << "_tpm";
	w->SetFileName(ss.str().c_str());
	w->Update();

	return EXIT_SUCCESS;
}
