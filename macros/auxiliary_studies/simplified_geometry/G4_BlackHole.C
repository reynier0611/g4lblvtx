#ifndef MACRO_G4_BLACKHOLE_C
#define MACRO_G4_BLACKHOLE_C

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4main/PHG4Reco.h>
#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

// ======================================================================================================================
void wrap_with_cylindrical_blackhole(PHG4Reco *g4Reco,double R,double zmin, double zmax, bool endcaps = true){
	/*
	   This code creates a cylinder of radius R and range in z from zmin to zmax
	   and makes these cylinders blackholes so that when a particle reaches them, these
	   surfaces absorb them and don't let these particles loop back into the geometry.
	   The values R, zmin, and zmax must be in cm.
	*/

	double thinness = 0.1;  // black hole thickness, needs to be taken into account for the z positions
	double hole = 0;	// If beampipe needs to be fed through a hole

	PHG4CylinderSubsystem * cyl_BH;

	cyl_BH = new PHG4CylinderSubsystem("BH_BARR",0);
	cyl_BH->set_double_param("length"   , zmax - zmin);
	cyl_BH->set_double_param("radius"   , R          );
	cyl_BH->set_double_param("thickness", thinness   );
	cyl_BH->set_string_param("material" , "G4_AIR"   );
	cyl_BH->SuperDetector("BH_BARR");
	cyl_BH->BlackHole();
	cyl_BH->set_color(0.9,0.5,1,0.5);
	g4Reco->registerSubsystem(cyl_BH);

	if(endcaps){
		cyl_BH = new PHG4CylinderSubsystem("BH_FOR",0);
		cyl_BH->set_double_param("length"   , thinness        );
		cyl_BH->set_double_param("radius"   , hole            );
		cyl_BH->set_double_param("thickness", R - hole        ); 
		cyl_BH->set_string_param("material" , "G4_AIR"        );
		cyl_BH->set_double_param("place_z"  , zmax+thinness/2.);
		cyl_BH->SuperDetector("BH_FOR");
		cyl_BH->BlackHole();
		cyl_BH->set_color(0.9,0.5,1,0.5);
		g4Reco->registerSubsystem(cyl_BH);

		cyl_BH = new PHG4CylinderSubsystem("BH_BACK",0);
		cyl_BH->set_double_param("length"   , thinness        );
		cyl_BH->set_double_param("radius"   , hole            );
		cyl_BH->set_double_param("thickness", R - hole        );
		cyl_BH->set_string_param("material" , "G4_AIR"        );
		cyl_BH->set_double_param("place_z"  , zmin-thinness/2.);
		cyl_BH->SuperDetector("BH_BACK");
		cyl_BH->BlackHole();
		cyl_BH->set_color(0.9,0.5,1,0.5);
		g4Reco->registerSubsystem(cyl_BH);
	}
}
#endif
