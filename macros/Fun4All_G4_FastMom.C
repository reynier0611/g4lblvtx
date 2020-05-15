#pragma once
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
#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>

#include <g4lblvtx/G4LBLVtxSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>
#include <g4lblvtx/TrackFastSimEval.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_FastMom(
			int nEvents = -1,
			const char *outputFile = "out_allSi",
			const char *genpar = "pi-")
{
	bool use_particle_gen = true;
        bool use_particle_gun = false;

// projections
	string projname = "DIRC";
        double projradius = 85.;

        if( (use_particle_gen&&use_particle_gun) && (!use_particle_gen&&!use_particle_gun) ){ cout << "Set one and only one variable above to true" << endl; exit(0);}
	
	cout << "Particle that will be generated: " << std::string(genpar) << endl;

	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// If you want to fix the random seed for reproducibility
	// recoConsts *rc = recoConsts::instance();
	// rc->set_IntFlag("RANDOMSEED", 12345);
	
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string(genpar));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);			// Vertex generation range
	gen->set_mom_range(1,10.);	// Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(-1,1);		// Detector coverage corresponds to |η|< 4
	gen->set_phi_range(0.,2.*TMath::Pi());
	// --------------------------------------------------------------------------------------
	// Particle Gun Setup
	PHG4ParticleGun *gun = new PHG4ParticleGun();
	gun->set_name(std::string(genpar));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ...
	gun->set_vtx(0,0,0);
	gun->set_mom(0,1,0);
	// --------------------------------------------------------------------------------------
	if     (use_particle_gen){se->registerSubsystem(gen); cout << "Using particle generator" << endl;}
	else if(use_particle_gun){se->registerSubsystem(gun); cout << "Using particle gun"       << endl;}
	
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	//g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
	//g4Reco->set_field_rescale(-1.4/1.5);
	float B_T = 1.5; // Magnetic Field [T]
	g4Reco->set_field(B_T);
	//g4Reco->SetPhysicsList("FTFP_BERT_HP");
	// ======================================================================================================
	// Loading All-Si Tracker from dgml file
	AllSiliconTrackerSubsystem *allsili = new AllSiliconTrackerSubsystem();
	//allsili->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/AllSiliconTracker/FAIRGeom.gdml");
	allsili->set_string_param("GDMPath","FAIRGeom.gdml"); 

	allsili->AddAssemblyVolume("VST");	// Barrel
	allsili->AddAssemblyVolume("FST");	// Forward disks
	allsili->AddAssemblyVolume("BST");	// Backward disks
	allsili->AddAssemblyVolume("BEAMPIPE");	// Beampipe

	// this is for plotting single logical volumes for debugging and geantino scanning they end up at the center, you can plot multiple
	// but they end up on top of each other. They cannot coexist with the assembly volumes, the code will quit if you try to use both.
	// allsili->AddLogicalVolume("BstContainerVolume04");
	// allsili->AddLogicalVolume("FstContainerVolume00");
	// allsili->AddLogicalVolume("FstChipAssembly37");
	// allsili->AddLogicalVolume("VstStave00");
	
	allsili->SuperDetector("LBLVTX");
	allsili->SetActive();          // this saves hits in the MimosaCore volumes
	allsili->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(allsili);

	// ======================================================================================================

	PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem(projname,0);
	cyl->set_double_param("length", 400);
	cyl->set_double_param("radius", projradius); // dirc radius
	cyl->set_double_param("thickness", 0.1); // needs some thickness
	cyl->set_string_param("material", "G4_AIR");
	cyl->SetActive(1);
	cyl->SuperDetector(projname);
	cyl->BlackHole();
        cyl->set_color(1,0,0,0.7); //reddish

	g4Reco->registerSubsystem(cyl);
	// ======================================================================================================

	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	// ======================================================================================================
	// Fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	char nodename[100];
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("SVTX");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	for (int i=10; i<16; i++){ // CENTRAL BARREL
		sprintf(nodename,"G4HIT_LBLVTX_CENTRAL_%d", i);
		kalman->add_phg4hits(
				nodename,				// const std::string& phg4hitsNames
				PHG4TrackFastSim::Cylinder,		// const DETECTOR_TYPE phg4dettype
				999.,					// radial-resolution [cm] (this number is not used in cylindrial geometry)
				5.8e-4,					// azimuthal (arc-length) resolution [cm]
				5.8e-4,					// longitudinal (z) resolution [cm]
				1,					// efficiency (fraction)
				0					// hit noise
				);
	}
	for (int i=20; i<25; i++){ // FORWARD DISKS
		sprintf(nodename,"G4HIT_LBLVTX_FORWARD_%d", i);
		kalman->add_phg4hits(
				nodename,                          	// const std::string& phg4hitsNames
				PHG4TrackFastSim::Vertical_Plane,  	// const DETECTOR_TYPE phg4dettype
				5.8e-4,                            	// radial-resolution [cm]
				5.8e-4,                            	// azimuthal (arc-length) resolution [cm]
				999.,                              	// longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
				1,                                 	// efficiency (fraction)
				0                                  	// hit noise
				);
	}
	for (int i=30; i<35; i++){ // BACKWARD DISKS
		sprintf(nodename,"G4HIT_LBLVTX_BACKWARD_%d", i);
		kalman->add_phg4hits(
				nodename,                          	// const std::string& phg4hitsNames
				PHG4TrackFastSim::Vertical_Plane,  	// const DETECTOR_TYPE phg4dettype
				5.8e-4,                            	// radial-resolution [cm]
				5.8e-4,                            	// azimuthal (arc-length) resolution [cm]
				999.,                              	// longitudinal (z) resolution [cm] (this number is not used in vertical plane geometry)
				1,                                 	// efficiency (fraction)
				0                                  	// hit noise
				);
	}
// projection on cylinder with 80cm radius
	kalman->add_cylinder_state(projname, projradius);
// projection on vertical plane at z=40cm
	kalman->add_zplane_state("MYZPLANE", 40.);
	se->registerSubsystem(kalman);
	// -----------------------------------------------------
	// INFO: The resolution numbers above correspond to:
	// 20e-4/sqrt(12) cm = 5.8e-4 cm, to simulate 20x20 um

	// ======================================================================================================
	TrackFastSimEval *fast_sim_eval = new TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(TString(outputFile)+Form("_B_%.1fT",B_T)+"_FastTrackingEval.root");
	fast_sim_eval->AddProjection(projname);
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	SimpleNtuple *hits = new SimpleNtuple("Hits");
	//  hits->AddNode("ABSORBER_LBLVTX",0); // hits in the passive volumes
	for (int i = 10; i < 16; i++){sprintf(nodename, "LBLVTX_CENTRAL_%d", i);	hits->AddNode(nodename, i);} // hits in the  MimosaCore volumes
	for (int i = 20; i < 25; i++){sprintf(nodename, "LBLVTX_FORWARD_%d", i);	hits->AddNode(nodename, i);} // hits in the  MimosaCore volumes
	for (int i = 30; i < 35; i++){sprintf(nodename, "LBLVTX_BACKWARD_%d",i);	hits->AddNode(nodename, i);} // hits in the  MimosaCore volumes
	se->registerSubsystem(hits);

	///////////////////////////////////////////
	// IOManagers...
	///////////////////////////////////////////
	const std::string dst_name = std::string(outputFile)+"_G4LBLVtx.root";
	//Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",TString(outputFile)+"_G4LBLVtx.root");
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
	out->Verbosity(10);
	se->registerOutputManager(out);

	Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);
	if (nEvents <= 0)
	{
		return;
	}
	se->run(nEvents);
	se->End();
	delete se;
	gSystem->Exit(0);
}
