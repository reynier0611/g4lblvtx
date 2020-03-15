#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
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

#include <g4lblvtx/G4LBLVtxSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
#endif


  void Fun4All_G4_gdml(
    int nEvents = -1
    )
{

  gSystem->Load("libfun4all");
  gSystem->Load("libg4detectors");
  gSystem->Load("libg4testbench");
  gSystem->Load("libg4histos");

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
  gen->set_name("e-");
  gen->set_vtx(0, -1, 0);
  //  gen->set_eta_range(0.297, 0.303);
  gen->set_mom_range(1.0, 10.0);
  gen->set_z_range(0.,0.);
// this eta/phi range roughly covers the stave 
// the corners are still outside this guarantees that the stave gets hit every time)
  gen->set_eta_range(-3.3, 3.3);
  gen->set_phi_range(60./180.*TMath::Pi(),120./180.*TMath::Pi());
  se->registerSubsystem(gen);

  PHG4ParticleGun *gun = new PHG4ParticleGun();
//  gun->set_name("pi-");
  gun->set_name("geantino");
  //gun->set_name("proton");
  gun->set_vtx(0, 0, 0);
  gun->set_mom(0, 1, 0);
  //se->registerSubsystem(gun);
  PHG4Reco* g4Reco = new PHG4Reco();
  g4Reco->set_field(0);
  g4Reco->save_DST_geometry(false);
  //g4Reco->SetPhysicsList("FTFP_BERT_HP");

  G4LBLVtxSubsystem *gdml = new G4LBLVtxSubsystem("LBLVtx");
  // gdml->set_double_param("place_x",30.);
  // gdml->set_double_param("place_y",30.);
  gdml->set_string_param("GDMPath","mvtx_stave_v02.gdml");
  gdml->set_string_param("TopVolName","MVTXStave_StaveStruct");
  gdml->SetActive();
  gdml->SetAbsorberActive();
  gdml->SuperDetector("LBLVTX");
  g4Reco->registerSubsystem(gdml);

  // PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  // g4Reco->registerSubsystem(truth);

  se->registerSubsystem( g4Reco );

  SimpleNtuple *hits = new SimpleNtuple("Hits");
  hits->AddNode("LBLVTX",0);
  hits->AddNode("ABSORBER_LBLVTX",1);
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
  gdml->Print();
  se->End();
  //   std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);

}

