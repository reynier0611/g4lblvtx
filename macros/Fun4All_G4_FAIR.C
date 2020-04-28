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

#include <g4lblvtx/G4LBLVtxSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libeicdetectors.so)


void Fun4All_G4_FAIR(
  int nEvents = -1
  )
{

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
//  gSystem->Load("libnodetest");
  // IntegralReco *intreg = new IntegralReco();
  // se->registerSubsystem(intreg);
  //  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RANDOMSEED",12345);
  PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
  gen->set_name("geantino");
  gen->set_vtx(0, 0, 0);
  //  gen->set_eta_range(0.297, 0.303);
  gen->set_mom_range(1.0, 1.0);
  gen->set_z_range(0.,0.);
  gen->set_eta_range(-3.3, 3.3); // not sure about the eta range
  gen->set_phi_range(0./180.*TMath::Pi(),360./180.*TMath::Pi());
  se->registerSubsystem(gen);

  PHG4ParticleGun *gun = new PHG4ParticleGun();
//  gun->set_name("pi-");
  gun->set_name("geantino");
  //gun->set_name("proton");
  gun->set_vtx(0, 0, 0);
  gun->set_mom(0, 1, 0);
  //se->registerSubsystem(gun);
  PHG4Reco* g4Reco = new PHG4Reco();
  g4Reco->set_field(1.5);
  g4Reco->save_DST_geometry(false);
  //g4Reco->SetPhysicsList("FTFP_BERT_HP");

  AllSiliconTrackerSubsystem *allsili = new AllSiliconTrackerSubsystem();
  allsili->set_string_param("GDMPath","FAIRGeom.gdml");


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

  allsili->SetActive(); // this saves hits in the MimosaCore volumes
  allsili->SetAbsorberActive(); // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(allsili);

  // PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  // g4Reco->registerSubsystem(truth);

  se->registerSubsystem( g4Reco );

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

  // char outnode[100];
  // sprintf(outnode,"G4RootHit_%s",detector);
  // out->AddNode(outnode);

  Fun4AllInputManager *in = new Fun4AllDummyInputManager( "JADE");
  se->registerInputManager( in );
  if (nEvents <= 0)
  {
    return;
  }
  se->run(nEvents);
//  gdml->Print();
  se->End();
  //   std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);

}
