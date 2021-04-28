#ifndef MACRO_G4GEMCYL_C
#define MACRO_G4GEMCYL_C

#include "GlobalVariables.C"
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4main/PHG4Reco.h>
#include <string>

/*
macro written by R. Cruz-Torres (reynier@lbl.gov)
adapted from the existing plane GEMs in Fun4All
*/

R__LOAD_LIBRARY(libg4detectors.so)
// ======================================================================================================================
int make_GEM_barrel_section(string name, PHG4Reco *g4Reco, double Rpos, double zpos , double X_X0_perc = 0.652604){
	if(X_X0_perc < 0.652604){
		X_X0_perc = 0.652604;
		std::cout << "\nX/X0 must be >= 0.652604. Defaulting to 0.652604!\n\n";
	}

	PHG4CylinderSubsystem *cyl;
        const double mm = 0.1;
        const double um = 1e-4;
	// ---------------------------------------------------------------------------------------------------
	// Mini-TPC
	TString mat_tpc[] = {"G4_MYLAR","G4_METHANE","G4_GRAPHITE"};
	const int n_layers_tpc = sizeof(mat_tpc)/sizeof(*mat_tpc);
        double thick_tpc[] = {  25*um  ,      2     ,   10*um     };
        double cl1_tpc[] = {1.0,0.1,0.1};
        double cl2_tpc[] = {0.6,0.3,0.7};
        double cl3_tpc[] = {0.1,0.9,0.2};

	double last_r = Rpos;
	
	for (int ilayer = 0; ilayer < n_layers_tpc; ilayer++){
                cyl = new PHG4CylinderSubsystem(name+"_miniTPC", ilayer);
                cyl->set_string_param("material" , (string)(mat_tpc[ilayer]));
                cyl->set_double_param("radius"   , last_r                   );
                cyl->set_double_param("thickness", thick_tpc[ilayer]        );
                cyl->set_double_param("place_z"  , 0                        );
                cyl->set_double_param("length"   , zpos                     );
		if(mat_tpc[ilayer]=="G4_METHANE")
                	cyl->SetActive();
                cyl->SuperDetector(name);
                cyl->set_color(cl1_tpc[ilayer],cl2_tpc[ilayer],cl3_tpc[ilayer]);
                g4Reco->registerSubsystem(cyl);
        
		last_r += thick_tpc[ilayer];
	}
	// ---------------------------------------------------------------------------------------------------
	// HBD
	// For the Hadron Blind Detector (HBD) details
        // see line 411 here: https://sphenix-collaboration.github.io/doxygen/d4/d18/PHG4SectorConstructor_8cc_source.html
        // GEM Copper 1.43 0.0005x6 64 0.134
        // GEM Kapton 28.6 0.005x3 64 0.034
        // GEM Copper 1.43 0.0005x6 64 0.134
        // GEM frames FR4 17.1 0.15x4 6.5 0.228
	const int nGEMLayers = 3; // for HBD
        TString mat_hbd1[]        = {"G4_Cu","G4_KAPTON","G4_Cu","G10"};
        const int n_layers_hbd1 = sizeof(mat_hbd1)/sizeof(*mat_hbd1);
        double thick_hbd1[]       = {0.0005 ,0.005      ,0.0005 ,0.15 };
	double perc_filled_hbd1[] = { .64   , .64       , .64   ,.065 };
        double cl1_hbd1[] = {1.0,0.0,1.0,0.9};
        double cl2_hbd1[] = {0.1,0.0,0.1,1.0};
        double cl3_hbd1[] = {0.1,0.5,0.1,0.8};

	for(int ngem = 0 ; ngem<nGEMLayers ; ngem++){
		for (int ilayer = 0; ilayer < n_layers_hbd1; ilayer++){
			cyl = new PHG4CylinderSubsystem(name, 10*(ngem+1)+ilayer);
                	cyl->set_string_param("material" , (string)(mat_hbd1[ilayer])                  );
                	cyl->set_double_param("radius"   , last_r                                      );
                	cyl->set_double_param("thickness", thick_hbd1[ilayer]*perc_filled_hbd1[ilayer] );
                	cyl->set_double_param("place_z"  , 0                                           );
                	cyl->set_double_param("length"   , zpos                                        ); 
                	cyl->SuperDetector(name);
                	cyl->set_color(cl1_hbd1[ilayer],cl2_hbd1[ilayer],cl3_hbd1[ilayer]);
                	g4Reco->registerSubsystem(cyl);
              
			last_r += thick_hbd1[ilayer];
		}
	}
	// -------------------------------------
	// PCB Kapton 28.6 0.005 100 0.017
	// PCB Copper 1.43 0.0005 80 0.028
	// Facesheet FR4 17.1 0.025x2 100 0.292
	TString mat_hbd2[] = {"G4_KAPTON","G4_Cu","G10"};
	const int n_layers_hbd2 = sizeof(mat_hbd2)/sizeof(*mat_hbd2);
        double thick_hbd2[]       = {0.005      ,0.0005 ,0.025*2};
        double perc_filled_hbd2[] = { 1         , .80   , 1     };
        double cl1_hbd2[] = {0.0,1.0,0.9};
        double cl2_hbd2[] = {0.0,0.1,1.0};
        double cl3_hbd2[] = {0.5,0.1,0.8};

	thick_hbd2[1] = thick_hbd2[1] + 1.43*(X_X0_perc-0.652604)/100./perc_filled_hbd2[1];

        for (int ilayer = 0; ilayer < n_layers_hbd2; ilayer++){
                cyl = new PHG4CylinderSubsystem(name+"_hbd", 40+ilayer);
                cyl->set_string_param("material" , (string)(mat_hbd2[ilayer])                  );
                cyl->set_double_param("radius"   , last_r                                      );
                cyl->set_double_param("thickness", thick_hbd2[ilayer]*perc_filled_hbd2[ilayer] );
                cyl->set_double_param("place_z"  , 0                                           );
                cyl->set_double_param("length"   , zpos                                        );
                cyl->SuperDetector(name);
                cyl->set_color(cl1_hbd2[ilayer],cl2_hbd2[ilayer],cl3_hbd2[ilayer]);
                g4Reco->registerSubsystem(cyl);

                last_r += thick_hbd2[ilayer];
        }
	// ---------------------------------------------------------------------------------------------------
	return 0;
}
// ======================================================================================================================
#endif
