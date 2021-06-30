/*
   ================================================================================================================
   The purpose of this code is to have a version of the all-silicon tracker that is simplified so that we can study
   different variations of the geometry quickly. Specifically, I wrote this code to study the impact that changing
   the material budget of different regions of the detector would have on different resolutions.
   ================================================================================================================
   */
#pragma once
#include <phgenfit/Track.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4histos/G4HitNtuple.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>
#include <g4lblvtx/EicFRichSubsystem.h>			// Forward RICH
#include <g4lblvtx/PHG4ParticleGenerator_flat_pT.h>	// Flat-pT generator
#include <g4lblvtx/AllSi_Al_support_Subsystem.h>	// Aluminum cone
#include "G4_BlackHole.C"				// Blackhole
#include "G4_Pipe_EIC.C"				// Beampipe
#include "G4_GEM.C"					// GEMs
#include "G4_DIRC_SMALL.C"				// DIRC

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_simplified_v2(
			int nEvents        = -1    ,	// number of events
			bool include_RICH  = false ,	// if true, RICH material will be included
			double GEM_res     = -1.   ,	// um, if > 0 forward, backward GEMs will be included
			int nDircSectors   = -1    ,	// Number of Quartz bars in the DIRC
			int magnetic_field = 5     ,	// Magnetic field setting
			bool do_projections = false,	// Projections
			TString out_name   = "out" )	// output filename
{	
	// ======================================================================================================
	// Input from the user
	double vtx_matBud  = 0.05; // % X/X0 (material budget of vertexing layers)
        double barr_matBud = 0.55; // % X/X0 (material budget of middle layers)
        double disk_matBud = 0.24; // % X/X0 (material budget of disk layers)
	const int particle_gen = 1;// 1 = particle generator, 2 = particle gun, 3 = simple event generator, 4 = pythia8 e+p collision, 5 = particle generator flat in pT
	double pix_size_vtx = 10.; // um - size of pixels in vertexing layers
	double pix_size_bar = 10.; // um - size of pixels in barrel layers
	double pix_size_dis = 10.; // um - size of pixels in disk layers
	bool use_blackhole = false;
	const int nDisks_per_side = 5;
	double pmin =  0.; // GeV/c
	double pmax = 30.; // GeV/c
	double eta_min = -3.5;
	double eta_max =  3.5;
	// ======================================================================================================
	// Parameters for projections
	string projname1   = "DIRC";            // Cylindrical surface object name
	double projradius1 = 50.;               // [cm] 
	double length1     = 400.;              // [cm]
	// ---
	double thinness    = 0.1;               // black hole thickness, needs to be taken into account for the z positions
	// ---
	string projname2   = "FOR";             // Forward plane object name
	double projzpos2   = 130+thinness/2.;   // [cm]
	double projradius2 = 50.;               // [cm]
	// ---
	string projname3   = "BACK";            // Backward plane object name
	double projzpos3   = -(130+thinness/2.);// [cm]
	double projradius3 = 50.;               // [cm] 
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// If you want to fix the random seed for reproducibility
	// recoConsts *rc = recoConsts::instance();
	// rc->set_IntFlag("RANDOMSEED", 12345);
	// ======================================================================================================
	if(nDircSectors>0)
		DIRCInit();
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string("pi-"));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);			// Vertex generation range
	gen->set_mom_range(pmin,pmax);		// Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(eta_min,eta_max);
	gen->set_phi_range(0,2.*TMath::Pi());
	// --------------------------------------------------------------------------------------
	// Particle generator flat in pT
	PHG4ParticleGenerator_flat_pT *gen_pT = new PHG4ParticleGenerator_flat_pT();
	gen_pT->set_name(std::string("pi-"));     // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen_pT->set_vtx(0,0,0);                    // Vertex generation range
	gen_pT->set_pT_range(pmin,pmax);         // Momentum generation range in GeV/c
	gen_pT->set_z_range(0.,0.);
	gen_pT->set_eta_range(eta_min,eta_max);               // Detector coverage corresponds to |η|< 4
	gen_pT->set_phi_range(0.,2.*TMath::Pi());
	// ======================================================================================================
	if     (particle_gen==1){se->registerSubsystem(   gen); cout << "Using particle generator"             << endl;}
	else if(particle_gen==5){se->registerSubsystem(gen_pT); cout << "Using particle generator flat in pT"  << endl;}
	else{ cout << "Particle generator option requested has not been implemented. Bailing out!" << endl; exit(0); }
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	//g4Reco->SetWorldMaterial("G4_Galactic");	
	// ======================================================================================================
	// Magnetic field setting
	TString B_label;
	if(magnetic_field==1){          // uniform 1.5T
		B_label = "_B_1.5T";
		g4Reco->set_field(1.5);
	}
	else if(magnetic_field==2){     // uniform 3.0T
		B_label = "_B_3.0T";
		g4Reco->set_field(3.0);
	}
	else if(magnetic_field==3){     // sPHENIX 1.4T map
		B_label = "_B_BaBar";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
		g4Reco->set_field_rescale(-1.4/1.5);
	}
	else if(magnetic_field==4){     // Beast 3.0T map
		B_label = "_B_BeAST";
		g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);
	}
	else if(magnetic_field==5){
		B_label = "_B_ATHENA_210507";
		TString path_to_map = "/global/homes/r/reynier/Singularity/BeastMagneticField/data/EIC_Magnetic_Field_Map_2021_05_07_radial_coords_[cm]_[T].120000.lines.Bmap";
		g4Reco->set_field_map(string(path_to_map), PHFieldConfig::kFieldBeast);
	}
	else if(magnetic_field==6){
		B_label = "_B_ATHENA_210528";
		TString path_to_map = "/global/homes/r/reynier/Singularity/BeastMagneticField/data/EIC_v.0.1.0_Magnetic_Field_Map_2021_05_28_radial_coords_[cm]_[T].401301.line.Bmap";
		g4Reco->set_field_map(string(path_to_map), PHFieldConfig::kFieldBeast);
	}
	else{                           // The user did not provide a valid B field setting
		cout << "User did not provide a valid magnetic field setting. Set 'magnetic_field'. Bailing out!" << endl;
	}	
	// ======================================================================================================
	// Detector setup
	PHG4CylinderSubsystem *cyl;
	//---------------------------
	// Vertexing
	double si_vtx_r_pos[] = {3.30,5.70};
	const int nVtxLayers = sizeof(si_vtx_r_pos)/sizeof(*si_vtx_r_pos);
	double si_z_vtxlength[] = {30.,30.};
	double si_thick_vtx = vtx_matBud/100.*9.37;

	for (int ilayer = 0; ilayer < nVtxLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
		cyl->set_string_param("material" , "G4_Si"               );
		cyl->set_double_param("radius"   , si_vtx_r_pos[ilayer]  );
		cyl->set_double_param("thickness", si_thick_vtx          );
		cyl->set_double_param("place_z"  , 0                     );
		cyl->set_double_param("length"   , si_z_vtxlength[ilayer]);
		cyl->SetActive();
		cyl->SuperDetector("SVTX");
		cyl->set_color(0,0.8,0.1);
		g4Reco->registerSubsystem(cyl);
	}
	//---------------------------
	// Barrel
	double si_r_pos[] = {21.,22.68,39.3,43.23};
	const int nTrckLayers = sizeof(si_r_pos)/sizeof(*si_r_pos);
	double si_z_length[] = {54.,60.,105.,114.};
	double si_thick_bar = barr_matBud/100.*9.37;

	for (int ilayer = 0; ilayer < nTrckLayers ; ilayer++){
		cyl = new PHG4CylinderSubsystem("BARR", ilayer);
		cyl->set_string_param("material" , "G4_Si"            );
		cyl->set_double_param("radius"   , si_r_pos[ilayer]   );
		cyl->set_double_param("thickness", si_thick_bar       );
		cyl->set_double_param("place_z"  , 0                  );
		cyl->set_double_param("length"   , si_z_length[ilayer]);
		cyl->SetActive();
		cyl->SuperDetector("BARR");
		cyl->set_color(0,0.5,1);
		g4Reco->registerSubsystem(cyl);	
	}
	//---------------------------
	// Disks
	const double z_min_d = 25.; // cm (z of the first disk)
	const double z_max_d = 121.; // cm (z of the last disk)
	const double disk_to_disk_distance = (z_max_d-z_min_d)/((float)nDisks_per_side-1.);
	const int nDisks = nDisks_per_side*2;
	double si_z_pos[nDisks] = {0};
	double si_r_max[nDisks] = {0};
	double si_r_min[nDisks] = {0};
	double si_thick_disk = disk_matBud/100.*9.37;

	for(int i = 0        ; i < nDisks_per_side   ; i++) si_z_pos[i] = -z_max_d + (float)i*disk_to_disk_distance;	
	for(int i = nDisks-1 ; i > nDisks_per_side-1 ; i--) si_z_pos[i] = z_min_d + (float)(i-nDisks_per_side)*disk_to_disk_distance;

	for(int i = 0 ; i < nDisks ; i++){
		si_r_max[i] = TMath::Min(43.23,18.5*abs(si_z_pos[i])/si_z_pos[nDisks_per_side]);
		if(si_z_pos[i]>66.8&&si_z_pos[i]>0) si_r_min[i] = (0.05025461*si_z_pos[i]-0.180808);
		else if(si_z_pos[i]>0) si_r_min[i] = 3.18;
		else if(si_z_pos[i]<-79.8&&si_z_pos[i]<0) si_r_min[i] = (-0.0297039*si_z_pos[i]+0.8058281);
		else si_r_min[i] = 3.18;
		si_r_max[i] -= si_r_min[i];
	}

	for (int ilayer = 0; ilayer < nDisks ; ilayer++){
		cyl = new PHG4CylinderSubsystem("FBVS", ilayer);
		cyl->set_string_param("material" , "G4_Si"         );
		cyl->set_double_param("radius"   , si_r_min[ilayer]);
		cyl->set_double_param("thickness", si_r_max[ilayer]);
		cyl->set_double_param("place_z"  , si_z_pos[ilayer]);
		cyl->set_double_param("length"   , si_thick_disk   );
		cyl->SetActive();
		cyl->SuperDetector("FBST");
		cyl->set_color(1,0,0);
		g4Reco->registerSubsystem(cyl);
	}
	//---------------------------
	// Black hole to suck loopers out of their misery
	double BH_r = si_r_pos[nTrckLayers-1]+2;
	double BH_zmin = si_z_pos[0]-2;
	double BH_zmax = si_z_pos[sizeof(si_z_pos)/sizeof(*si_z_pos)-1]+2;
	if(use_blackhole)
		wrap_with_cylindrical_blackhole(g4Reco,BH_r,BH_zmin,BH_zmax);
	//---------------------------
	// Beam pipe
	PipeInit(); // Load beampipe from Fun4All rather than from gdml file
	double pipe_radius = 0;
	pipe_radius = Pipe(g4Reco,pipe_radius,true);
	//---------------------------
	// Al Support Structure
	AllSi_Al_support_Subsystem *Al_supp = new AllSi_Al_support_Subsystem("Al_supp");
	g4Reco->registerSubsystem(Al_supp);
	//---------------------------
	// DIRC
	if(nDircSectors>0)
		double dirc_out_skin = DIRCSetup(g4Reco,nDircSectors);	
	//---------------------------
	// forward RICH
	if(include_RICH){
		EicFRichSubsystem *RICH = new EicFRichSubsystem("RICH");
		g4Reco->registerSubsystem(RICH);
	}
	//---------------------------
	// Forward - Backward GEMs
	if(GEM_res>0){
		EGEM_Init();                            // Loading backward GEM geometry
		FGEM_Init();                            // Loading forward GEM geometry
		EGEMSetup(g4Reco);
		FGEMSetup(g4Reco);
	}
	//---------------------------
	if(do_projections){
		PHG4CylinderSubsystem *cyl;
		cyl = new PHG4CylinderSubsystem(projname1,0);
		cyl->set_double_param("length", length1);
		cyl->set_double_param("radius", projradius1); // dirc radius
		cyl->set_double_param("thickness", 0.1); // needs some thickness
		cyl->set_string_param("material", "G4_AIR");
		cyl->SetActive(1);
		cyl->SuperDetector(projname1);
		cyl->BlackHole();
		cyl->set_color(1,0,0,0.7); //reddish
		g4Reco->registerSubsystem(cyl);

		cyl = new PHG4CylinderSubsystem(projname2,0);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 2); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius2); // 
		cyl->set_string_param("material", "G4_AIR");
		cyl->set_double_param("place_z", projzpos2);
		cyl->SetActive(1);
		cyl->SuperDetector(projname2);
		cyl->BlackHole();
		cyl->set_color(0,1,1,0.3); //reddish
		g4Reco->registerSubsystem(cyl);

		cyl = new PHG4CylinderSubsystem(projname3,0);
		cyl->set_double_param("length", thinness);
		cyl->set_double_param("radius", 2); // beampipe needs to fit here
		cyl->set_double_param("thickness", projradius3); // 
		cyl->set_string_param("material", "G4_AIR");
		cyl->set_double_param("place_z", projzpos3);
		cyl->SetActive(1);
		cyl->SuperDetector(projname3);
		cyl->BlackHole();
		cyl->set_color(0,1,1,0.3); //reddish
		g4Reco->registerSubsystem(cyl);
	}
	//---------------------------
	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	// ======================================================================================================
	// fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//---------------------------
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("BARR");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	// add Vertexing Layers
	kalman->add_phg4hits(
			"G4HIT_SVTX",				// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,					// radial-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),		// azimuthal-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),		// z-resolution [cm]
			1,					// efficiency,
			0					// noise hits
			);

	// add Barrel Layers
	kalman->add_phg4hits(
			"G4HIT_BARR",                   	// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,                           	// radial-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// azimuthal-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      	// z-resolution [cm]
			1,                              	// efficiency,
			0                               	// noise hits
			);

	// add Disk Layers
	kalman->add_phg4hits(
			"G4HIT_FBST",				// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Vertical_Plane,
			pix_size_dis/10000./sqrt(12.),		// radial-resolution [cm]
			pix_size_dis/10000./sqrt(12.),		// azimuthal-resolution [cm]
			999.,                       		// z-resolution [cm]
			1,                          		// efficiency,
			0                           		// noise hits
			);	
	// add forward and backward GEMs
	if(GEM_res>0){
		// BACKWARD GEM
		kalman->add_phg4hits("G4HIT_EGEM",              // const std::string& phg4hitsNames,
				PHG4TrackFastSim::Vertical_Plane,               // const DETECTOR_TYPE phg4dettype,
				GEM_res/10000.,                                 // const float radres,
				GEM_res/10000.,                                 // const float phires,
				999.,                                           // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
				1,                                              // const float eff,
				0                                               // const float noise
				);
		// FORWARD GEM2
		kalman->add_phg4hits("G4HIT_FGEM",              // const std::string& phg4hitsNames,
				PHG4TrackFastSim::Vertical_Plane,               // const DETECTOR_TYPE phg4dettype,
				GEM_res/10000.,                                 // const float radres,
				GEM_res/10000.,                                 // const float phires,
				999.,                                           // longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
				1,                                              // const float eff,
				0                                               // const float noise
				);	
	}

	// Projections  
	if(do_projections){
		kalman->add_cylinder_state(projname1, projradius1);     // projection on cylinder (DIRC)
		kalman->add_zplane_state  (projname2, projzpos2  );     // projection on vertical planes
		kalman->add_zplane_state  (projname3, projzpos3  );     // projection on vertical planes
	}

	//kalman->Verbosity(10);
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_vertex_xy_resolution(0);
	kalman->set_vertex_z_resolution(0);
	kalman->enable_vertexing(false); // this is false by default
	kalman->set_vertex_min_ndf(2);

	se->registerSubsystem(kalman);

	// ======================================================================================================
	TString label_mat = Form("_AllSi_vbd_%.2f_%.2f_%.2f",vtx_matBud,barr_matBud,disk_matBud);
	TString label_RICH = "";
	TString label_GEM  = "";
	TString label_DIRC = "";

	if(include_RICH)
		label_RICH = "_RICH";

	if(GEM_res>0)
		label_GEM = Form("_GEM_res_%.1fum",GEM_res);

	if(nDircSectors>0)
		label_DIRC = Form("_DIRC_%i_sect",nDircSectors);

	std::string outputFile = (std::string)(out_name)+std::string(label_mat)+std::string(label_RICH)+std::string(label_GEM)+std::string(label_DIRC)+std::string(B_label)+"_FastSimEval.root";
	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(outputFile);
	if(do_projections){
                fast_sim_eval->AddProjection(projname1);
                fast_sim_eval->AddProjection(projname2);
                fast_sim_eval->AddProjection(projname3);
        }
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	// IOManagers...
	const std::string dst_name = std::string(out_name)+std::string(B_label)+"_G4LBLVtx.root";	
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
	out->Verbosity(0);
	se->registerOutputManager(out);

	Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);

	if (nEvents <= 0) return;

	se->run(nEvents);
	se->End();
	delete se;

	gSystem->Exit(0);
}
