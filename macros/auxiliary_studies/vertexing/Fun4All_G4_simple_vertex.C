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

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4lblvtx.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void Fun4All_G4_simple_vertex(
			int nEvents = -1,		// number of events
			bool vtx_lyr_1 = true,
			bool vtx_lyr_2 = true,
			bool vtx_lyr_3 = true,
			double pmin = 0., // GeV/c
			double pmax = 30., // GeV/c
			double Bfield = 3.0, //T
			TString out_name = "out_vtx_study")	// output filename
{	
	TString outputFile = out_name+Form("_B_%.1fT",Bfield)+"_FastSimEval.root";
	double vtx_matBud = 0.05; //% X/X0
        double pix_size_vtx = 10.; // um - size of pixels in vertexing layers
        double pix_size_bar = 10.; // um - size of pixels in barrel layers
        double pix_size_dis = 10.; // um - size of pixels in disk layers
	// ======================================================================================================
	// Make the Server
	Fun4AllServer *se = Fun4AllServer::instance();
	// If you want to fix the random seed for reproducibility
	// recoConsts *rc = recoConsts::instance();
	// rc->set_IntFlag("RANDOMSEED", 12345);
	// ======================================================================================================
	// Particle Generator Setup
	PHG4ParticleGenerator *gen = new PHG4ParticleGenerator();
	gen->set_name(std::string("pi-"));	// geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... (currently passed as an input)
	gen->set_vtx(0,0,0);			// Vertex generation range
	gen->set_mom_range(pmin,pmax);		// Momentum generation range in GeV/c
	gen->set_z_range(0.,0.);
	gen->set_eta_range(1.,3.);
	gen->set_phi_range(0,2.*TMath::Pi());
	se->registerSubsystem(gen);
	// ======================================================================================================
	PHG4Reco *g4Reco = new PHG4Reco();
	g4Reco->SetWorldMaterial("G4_Galactic");	
	g4Reco->set_field(Bfield);
	// ======================================================================================================
	// Detector setup
	PHG4CylinderSubsystem *cyl;
	//---------------------------
	// Vertexing
	double si_vtx_r_pos[] = {3.64,4.45,5.26};
	const int nVtxLayers = sizeof(si_vtx_r_pos)/sizeof(*si_vtx_r_pos);
	double si_z_vtxlength[] = {14.,14.,14.};
	for(int i = 0 ; i < nVtxLayers ; i++) si_z_vtxlength[i] *= 3.;
	double si_thick_vtx = vtx_matBud/100.*9.37;

	for (int ilayer = 0; ilayer < nVtxLayers ; ilayer++){
		if(
				(ilayer==0&&vtx_lyr_1)||
				(ilayer==1&&vtx_lyr_2)||
				(ilayer==2&&vtx_lyr_3)
		  ){
			cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
			cyl->set_string_param("material" , "G4_Si"               );
			cyl->set_double_param("radius"   , si_vtx_r_pos[ilayer]  );
			cyl->set_double_param("thickness", si_thick_vtx          );
			cyl->set_double_param("place_z"  , 0                     );
			cyl->set_double_param("length"   , si_z_vtxlength[ilayer]);
			cyl->SetActive();
			cyl->SuperDetector("SVTX");
			g4Reco->registerSubsystem(cyl);
		}
	}
	//---------------------------
	// Barrel
	double si_r_pos[] = {21.,22.68,39.3,43.23};
	const int nTrckLayers = sizeof(si_r_pos)/sizeof(*si_r_pos);
	double si_z_length[] = {18.,20.,35.,38.};
	for(int i = 0 ; i < nTrckLayers ; i++) si_z_length[i] *= 3.;
	double si_thick_bar = 0.55/100.*9.37;

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
	double si_z_pos[] = {-121.,-97.,-73.,-49.,-25.,25.,49.,73.,97.,121.};
	double si_r_max[10] = {0};
	double si_r_min[10] = {0};
	double si_thick_disk = 0.3/100.*9.37;
	for(int i = 0 ; i < 10 ; i++){
		si_r_max[i] = TMath::Min(43.23,18.5*abs(si_z_pos[i])/si_z_pos[5]);

		if(si_z_pos[i]>66.8&&si_z_pos[i]>0) si_r_min[i] = (0.05025461*si_z_pos[i]-0.180808);
		else if(si_z_pos[i]>0) si_r_min[i] = 3.18;
		else if(si_z_pos[i]<-79.8&&si_z_pos[i]<0) si_r_min[i] = (-0.0297039*si_z_pos[i]+0.8058281);
		else si_r_min[i] = 3.18;
	}

	for (int ilayer = 0; ilayer < 10; ilayer++){
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
	// mid-rapidity beryllium pipe
	double be_pipe_radius = 3.1000;
	double be_pipe_thickness = 3.1762 - be_pipe_radius;  // 760 um for sPHENIX
	double be_pipe_length_plus = 66.8;                   // +z beam pipe extend.
	double be_pipe_length_neg = -79.8;                   // -z beam pipe extend.
	double be_pipe_length = be_pipe_length_plus - be_pipe_length_neg;
	double be_pipe_center = 0.5 * (be_pipe_length_plus + be_pipe_length_neg);

	cyl = new PHG4CylinderSubsystem("BE_PIPE", 1);
	cyl->set_double_param("radius", be_pipe_radius);
	cyl->set_int_param("lengthviarapidity", 0);
	cyl->set_double_param("length", be_pipe_length);
	cyl->set_double_param("place_z", be_pipe_center);
	cyl->set_string_param("material", "G4_Be");
	cyl->set_double_param("thickness", be_pipe_thickness);
	cyl->SuperDetector("PIPE");
	g4Reco->registerSubsystem(cyl);
	//---------------------------

	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	//---------------------------
	// fast pattern recognition and full Kalman filter
	// output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//---------------------------
	PHG4TrackFastSim *kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_sub_top_node_name("BARR");
	kalman->set_trackmap_out_name("SvtxTrackMap");

	// add Vertexing Layers
	kalman->add_phg4hits(
			"G4HIT_SVTX",			// const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,				// radial-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),	// azimuthal-resolution [cm]
			pix_size_vtx/10000./sqrt(12.),	// z-resolution [cm]
			1,				// efficiency,
			0				// noise hits
			);
	
	// add Barrel Layers
	kalman->add_phg4hits(
			"G4HIT_BARR",                   // const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder,
			999.,                           // radial-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      // azimuthal-resolution [cm]
			pix_size_bar/10000./sqrt(12.),      // z-resolution [cm]
			1,                              // efficiency,
			0                               // noise hits
			);
	
	//  add Disk Layers
	kalman->add_phg4hits(
			"G4HIT_FBST",               // const std::string& phg4hitsNames,
			PHG4TrackFastSim::Vertical_Plane,
			pix_size_dis/10000./sqrt(12.),       // radial-resolution [cm]
			pix_size_dis/10000./sqrt(12.),       // azimuthal-resolution [cm]
			999.,                       // z-resolution [cm]
			1,                          // efficiency,
			0                           // noise hits
			);	
	
	//kalman->Verbosity(10);
	kalman->set_use_vertex_in_fitting(false);
	kalman->set_vertex_xy_resolution(0);
	kalman->set_vertex_z_resolution(0);
	kalman->enable_vertexing(false); // this is false by default
	kalman->set_vertex_min_ndf(2);

	se->registerSubsystem(kalman);

	PHG4TrackFastSimEval *fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_filename(outputFile);
	se->registerSubsystem(fast_sim_eval);

	// ======================================================================================================
	// IOManagers...
	const std::string dst_name = std::string(out_name)+"_G4LBLVtx.root";	
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
