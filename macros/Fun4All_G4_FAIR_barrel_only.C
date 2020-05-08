#pragma once
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4histos/G4HitNtuple.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <phool/recoConsts.h>
#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <g4lblvtx/G4LBLVtxSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>

#include <g4eval/SvtxEvaluator.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
//R__LOAD_LIBRARY(libeicdetectors.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_FAIR_barrel_only(
	int nEvents = -1,
	const char *outputFile = "barrel_only"
){

	bool use_particle_gen = true;
	bool use_particle_gun = false;

	if( (use_particle_gen&&use_particle_gun) && (!use_particle_gen&&!use_particle_gun) ){ cout << "Set one and only one variable above to true" << endl; exit(0);}

	///////////////////////////////////////////
	// Make the Server
	//////////////////////////////////////////
	Fun4AllServer *se = Fun4AllServer::instance();

	recoConsts *rc = recoConsts::instance();
	//rc->set_IntFlag("RANDOMSEED",12345);
	
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name("pi-"); // geantino, pi-, proton, ...
	gen->set_vtx(0, 0, 0);
	gen->set_mom_range(1.0, 1.0); // Momentum generation range in GeV
	gen->set_z_range(0.,0.);
	gen->set_eta_range(-3.3, 3.3); // not sure about the eta range
	gen->set_phi_range(0./180.*TMath::Pi(),360./180.*TMath::Pi());
	
	// --------------------------------------------------------------------------------------
	// Particle Gun Setup
	PHG4ParticleGun *gun = new PHG4ParticleGun();
        gun->set_name("pi-");
        gun->set_vtx(0, 0, 0);
        gun->set_mom(0, 1, 0);
       
	// --------------------------------------------------------------------------------------
	if     (use_particle_gen){se->registerSubsystem(gen); cout << "Using particle generator" << endl;}
	else if(use_particle_gun){se->registerSubsystem(gun); cout << "Using particle gun"       << endl;}

	// ======================================================================================================
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.5);
	//g4Reco->save_DST_geometry(false);
	//g4Reco->SetPhysicsList("FTFP_BERT_HP");

	// ======================================================================================================
	// Loading All-Si Tracker from dgml file
	AllSiliconTrackerSubsystem *allsili = new AllSiliconTrackerSubsystem();
	allsili->set_string_param("GDMPath","FAIRGeom.gdml");

	allsili->AddAssemblyVolume("VST");
	// The two lines below are commented out until we have the necessary code to get the forward and backward regions working
	// allsili->AddAssemblyVolume("FST");
	// allsili->AddAssemblyVolume("BST");
	allsili->AddAssemblyVolume("BEAMPIPE");

	// this is for plotting single logical volumes for debugging
	// and geantino scanning they end up at the center, you can plot multiple 
	// but they end up on top of each other. They cannot coexist with the
	// assembly volumes, the code will quit if you try to use both
	// allsili->AddLogicalVolume("BstContainerVolume04");
	// allsili->AddLogicalVolume("FstContainerVolume00");
	// allsili->AddLogicalVolume("FstChipAssembly37");
	// allsili->AddLogicalVolume("VstStave00");
	allsili->SuperDetector("LBLVTX");

	allsili->SetActive(); // this saves hits in the MimosaCore volumes
	allsili->SetAbsorberActive(); // this saves hits in all volumes (in the absorber node)
	g4Reco->registerSubsystem(allsili);

	// ======================================================================================================
	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem( g4Reco );

	// ======================================================================================================
	// fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//---------------------------
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("LBLVTX");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	kalman->add_phg4hits(
	"G4HIT_LBLVTX",             // const std::string& phg4hitsNames,
	PHG4TrackFastSim::Cylinder, // const DETECTOR_TYPE phg4dettype,
	300e-4,                     // radial-resolution [cm] (irrelevant for cylindrical geometry. Use any number)
	5e-4,                       // azimuthal-resolution [cm] (currently set to value from Alice)
	5e-4,                       // z-resolution [cm] (currently set to value from Alice)
	1,                          // efficiency,
	0                           // noise hits
	);
	se->registerSubsystem(kalman);

	// ======================================================================================================
	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
        fast_sim_eval->set_filename(TString(outputFile)+"_FastTrackingEval.root");
        se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	SimpleNtuple *hits = new SimpleNtuple("Hits");
	hits->AddNode("LBLVTX",0); // hits in the  MimosaCore volumes
	hits->AddNode("ABSORBER_LBLVTX",1); // hits in the passive volumes
	se->registerSubsystem(hits);

	///////////////////////////////////////////
	// IOManagers...
	///////////////////////////////////////////

	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT","G4LBLVtx.root");
	out->Verbosity(10);
	se->registerOutputManager(out);

	Fun4AllInputManager *in = new Fun4AllDummyInputManager( "JADE");
	se->registerInputManager( in );
	if (nEvents <= 0){return;}
	se->run(nEvents);

	se->End();

	std::cout << "All done" << std::endl;
	delete se;
	gSystem->Exit(0);
}
