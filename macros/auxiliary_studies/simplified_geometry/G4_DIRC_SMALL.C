#ifndef MACRO_G4DIRC_SMALL_C
#define MACRO_G4DIRC_SMALL_C

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <cmath>

R__LOAD_LIBRARY(libg4detectors.so)

	/*
	  Macro written by Reynier Cruz-Torres (reynier@lbl.gov)
	  based on G4_DIRC.C by Jin Huang (jhuang@bnl.gov)
	  */

	namespace Enable
{
	bool DIRC_SMALL = false;
	bool DIRC_SMALL_OVERLAPCHECK = false;
}  // namespace Enable

//namespace G4DIRC_SMALL
//{
	//double z_end = z_shift - length / 2.;
	//double outer_skin_radius = 58;
//}  // namespace G4DIRC_SMALL

void DIRCInit()
{
	//BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4DIRC_SMALL::outer_skin_radius);
	//BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4DIRC_SMALL::z_start);
	//BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4DIRC_SMALL::z_end);
}

//! Babar DIRC_SMALL (Without most of support structure)
//! Ref: I. Adam et al. The DIRC particle identification system for the BaBar experiment.
//! Nucl. Instrum. Meth., A538:281-357, 2005. doi:10.1016/j.nima.2004.08.129.
double DIRCSetup(PHG4Reco *g4Reco,int nSectors = 12)
{
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::DIRC_SMALL_OVERLAPCHECK;

	// Preserving the length of the BaBar quartz slates requires quantizing the radius
	double R12 = 83.65;
	double l_quartz = 2.*R12*TMath::Sin(TMath::Pi()/12.);
	double radiator_R = l_quartz/2./TMath::Sin(TMath::Pi()/((float)nSectors));

	double scale = 12./((float)nSectors);

	double length = 400;
	double z_shift = -75;
	double z_start = z_shift + length / 2.;
	double z_end = z_shift - length / 2.;
	//double outer_skin_radius = 89.25/scale;
	double outer_skin_radius = TMath::Sqrt(l_quartz*l_quartz/4+radiator_R*radiator_R) + 2.843672369*(12./((float)nSectors));

	PHG4SectorSubsystem *dirc;
	dirc = new PHG4SectorSubsystem("DIRC_SMALL");
	dirc->get_geometry().set_normal_polar_angle(M_PI / 2);
	dirc->get_geometry().set_normal_start(radiator_R * PHG4Sector::Sector_Geometry::Unit_cm());
	dirc->get_geometry().set_min_polar_angle(atan2(radiator_R, z_start));
	dirc->get_geometry().set_max_polar_angle(atan2(radiator_R, z_end));
	dirc->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	dirc->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	dirc->get_geometry().set_material("Quartz");
	dirc->get_geometry().set_N_Sector(nSectors);
	dirc->OverlapCheck(OverlapCheck);
	dirc->get_geometry().AddLayer("Radiator", "Quartz", 1.7 * PHG4Sector::Sector_Geometry::Unit_cm(), true);
	g4Reco->registerSubsystem(dirc);

	PHG4CylinderSubsystem *cyl;

	//  The cylinder skins provide most of the strength
	//  and stiffness of the CST. The thickness of the inner
	//  and outer skins is 1.27 and 0.76 mm, respectively

	// Inner skin:
	cyl = new PHG4CylinderSubsystem("DIRC_SMALL_CST_Inner_Skin", 10);
	cyl->set_double_param("radius", 81.71/scale);
	cyl->set_double_param("length", length);
	cyl->set_string_param("material", "G4_Al");
	cyl->set_double_param("thickness", 0.127);
	cyl->set_double_param("place_x", 0.);
	cyl->set_double_param("place_y", 0.);
	cyl->set_double_param("place_z", z_shift);
	cyl->SetActive(0);
	cyl->SuperDetector("DIRC_SMALL");
	cyl->OverlapCheck(OverlapCheck);

	g4Reco->registerSubsystem(cyl);

	// Outer skin:
	cyl = new PHG4CylinderSubsystem("DIRC_SMALL_CST_Outer_Skin", 11);
	cyl->set_double_param("radius", outer_skin_radius - 0.076);
	cyl->set_double_param("length", length);
	cyl->set_string_param("material", "G4_Al");
	cyl->set_double_param("thickness", 0.076);
	cyl->set_double_param("place_x", 0.);
	cyl->set_double_param("place_y", 0.);
	cyl->set_double_param("place_z", z_shift);
	cyl->SetActive(0);
	cyl->SuperDetector("DIRC_SMALL");
	cyl->OverlapCheck(OverlapCheck);

	g4Reco->registerSubsystem(cyl);
	
	// Done
	return outer_skin_radius;
}
#endif
