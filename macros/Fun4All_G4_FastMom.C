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

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_FastMom(
    int nEvents = -1)
{
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
// if you want to fix the random seed for reproducibility
//  recoConsts *rc = recoConsts::instance();
//  rc->set_IntFlag("RANDOMSEED", 12345);
  PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
  gen->set_name("pi-");
  gen->set_vtx(0, 0, 0);
  //  gen->set_eta_range(0.297, 0.303);
  gen->set_mom_range(1.0, 1.0);
  gen->set_z_range(0., 0.);
  gen->set_eta_range(-0.2, 0.2);  // not sure about the eta range
  gen->set_phi_range(0. / 180. * TMath::Pi(), 360. / 180. * TMath::Pi());
  se->registerSubsystem(gen);

  PHG4ParticleGun *gun = new PHG4ParticleGun();
  //  gun->set_name("pi-");
  gun->set_name("geantino");
  //gun->set_name("proton");
  gun->set_vtx(0, 0, 0);
  gun->set_mom(0, 1, 0);
  //se->registerSubsystem(gun);
  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root"), PHFieldConfig::kField2D);
  g4Reco->set_field_rescale(-1.4/1.5);
  //g4Reco->SetPhysicsList("FTFP_BERT_HP");

  AllSiliconTrackerSubsystem *allsili = new AllSiliconTrackerSubsystem();
  allsili->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/AllSiliconTracker/FAIRGeom.gdml");

  allsili->AddAssemblyVolume("VST");
  allsili->AddAssemblyVolume("FST");
  allsili->AddAssemblyVolume("BST");
  allsili->AddAssemblyVolume("BEAMPIPE");

  // this is for plotting single logical volumes for debugging
  // and geantino scanning they end up at the center, you can plot multiple
  // but they end up on top of each other. They cannot coexist with the
  // assembly volumes, the code will quit if you try to use both
  //    allsili->AddLogicalVolume("BstContainerVolume04");
  //allsili->AddLogicalVolume("FstContainerVolume00");
  //    allsili->AddLogicalVolume("FstChipAssembly37");
  //allsili->AddLogicalVolume("VstStave00");
  allsili->SuperDetector("LBLVTX");

  allsili->SetActive();          // this saves hits in the MimosaCore volumes
  allsili->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(allsili);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  se->registerSubsystem(g4Reco);
  char nodename[100];
  PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
  kalman->set_use_vertex_in_fitting(false);
  kalman->set_sub_top_node_name("SVTX");
  kalman->set_trackmap_out_name("SvtxTrackMap");
  for (int i=10; i<16; i++)
  {
    sprintf(nodename,"G4HIT_LBLVTX_CENTRAL_%d", i);
    kalman->add_phg4hits(
      nodename,                //      const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,
      300e-4,                      //       radial-resolution [cm]
      30e-4,                       //        azimuthal-resolution [cm]
      1,                           //      z-resolution [cm]
      1,                           //      efficiency,
      0                            //      noise hits
      );
  }
  // forward disk (+ 2 inner barrel layers)
  for (int i=20; i<25; i++)
  {
    sprintf(nodename,"G4HIT_LBLVTX_FORWARD_%d", i);
    kalman->add_phg4hits(
      nodename,                    //      const std::string& phg4hitsNames,
      PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
      1. / sqrt(12.),                    //      const float radres,
      70e-4,                             //      const float phires,
      100e-4,                            //      const float lonres,
      1,                                 //      const float eff,
      0                                  //      const float noise
      );
  }
  for (int i=30; i<35; i++)
  {
    sprintf(nodename,"G4HIT_LBLVTX_BACKWARD_%d", i);
    kalman->add_phg4hits(
      nodename,                    //      const std::string& phg4hitsNames,
      PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
      1. / sqrt(12.),                    //      const float radres,
      70e-4,                             //      const float phires,
      100e-4,                            //      const float lonres,
      1,                                 //      const float eff,
      0                                  //      const float noise
      );
  }

  se->registerSubsystem(kalman);

  PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
  fast_sim_eval->set_filename("FastTrackingEval.root");
  se->registerSubsystem(fast_sim_eval);

  SimpleNtuple *hits = new SimpleNtuple("Hits");
  //  hits->AddNode("ABSORBER_LBLVTX",0); // hits in the passive volumes
  for (int i = 10; i < 16; i++)
  {
    sprintf(nodename, "LBLVTX_CENTRAL_%d", i);
    hits->AddNode(nodename, i);  // hits in the  MimosaCore volumes
  }
  for (int i = 20; i < 25; i++)
  {
    sprintf(nodename, "LBLVTX_FORWARD_%d", i);
    hits->AddNode(nodename, i);  // hits in the  MimosaCore volumes
  }
  for (int i = 30; i < 35; i++)
  {
    sprintf(nodename, "LBLVTX_BACKWARD_%d", i);
    hits->AddNode(nodename, i);  // hits in the  MimosaCore volumes
  }
  se->registerSubsystem(hits);

  ///////////////////////////////////////////
  // IOManagers...
  ///////////////////////////////////////////

  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", "G4LBLVtx.root");
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
