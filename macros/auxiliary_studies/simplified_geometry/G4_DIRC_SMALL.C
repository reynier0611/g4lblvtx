#ifndef MACRO_G4DIRC_SMALL_C
#define MACRO_G4DIRC_SMALL_C

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>

#include <g4main/PHG4Reco.h>

#include <cmath>

R__LOAD_LIBRARY(libg4detectors.so)

	/*!
	 *  * \file G4_DIRC_SMALL.C
	 *   * \brief Macro setting up the barrel DIRC_SMALL
	 *    * \author Jin Huang <jhuang@bnl.gov>
	 *     * \version $Revision: 1.3 $
	 *      * \date $Date: 2013/10/09 01:08:17 $
	 *       */

	namespace Enable
{
	bool DIRC_SMALL = false;
	bool DIRC_SMALL_OVERLAPCHECK = false;
}  // namespace Enable

namespace G4DIRC_SMALL
{
	double radiator_R = 49.90;
	double length = 400;
	double z_shift = -75;  //115
	double z_start = z_shift + length / 2.;
	double z_end = z_shift - length / 2.;
	double outer_skin_radius = 58;
}  // namespace G4DIRC_SMALL

void DIRCInit()
{
	BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4DIRC_SMALL::outer_skin_radius);
	BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4DIRC_SMALL::z_start);
	BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, G4DIRC_SMALL::z_end);
}

//! Babar DIRC_SMALL (Without most of support structure)
//! Ref: I. Adam et al. The DIRC particle identification system for the BaBar experiment.
//! Nucl. Instrum. Meth., A538:281-357, 2005. doi:10.1016/j.nima.2004.08.129.
double DIRCSetup(PHG4Reco *g4Reco)
{
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::DIRC_SMALL_OVERLAPCHECK;

	PHG4SectorSubsystem *dirc;
	dirc = new PHG4SectorSubsystem("DIRC_SMALL");
	dirc->get_geometry().set_normal_polar_angle(M_PI / 2);
	dirc->get_geometry().set_normal_start(49.90 * PHG4Sector::Sector_Geometry::Unit_cm());
	dirc->get_geometry().set_min_polar_angle(atan2(G4DIRC_SMALL::radiator_R, G4DIRC_SMALL::z_start));
	dirc->get_geometry().set_max_polar_angle(atan2(G4DIRC_SMALL::radiator_R, G4DIRC_SMALL::z_end));
	dirc->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	dirc->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::FlatEdge());
	dirc->get_geometry().set_material("Quartz");
	dirc->get_geometry().set_N_Sector(7);
	dirc->OverlapCheck(OverlapCheck);
	dirc->get_geometry().AddLayer("Radiator", "Quartz", 1.7 * PHG4Sector::Sector_Geometry::Unit_cm(), true);
	g4Reco->registerSubsystem(dirc);

	PHG4CylinderSubsystem *cyl;

	//  The cylinder skins provide most of the strength
	//  and stiffness of the CST. The thickness of the inner
	//  and outer skins is 1.27 and 0.76 mm, respectively

	// Inner skin:
	cyl = new PHG4CylinderSubsystem("DIRC_SMALL_CST_Inner_Skin", 10);
	cyl->set_double_param("radius", 47.96);
	cyl->set_double_param("length", G4DIRC_SMALL::length);
	cyl->set_string_param("material", "G4_Al");
	cyl->set_double_param("thickness", 0.127);
	cyl->set_double_param("place_x", 0.);
	cyl->set_double_param("place_y", 0.);
	cyl->set_double_param("place_z", G4DIRC_SMALL::z_shift);
	cyl->SetActive(0);
	cyl->SuperDetector("DIRC_SMALL");
	cyl->OverlapCheck(OverlapCheck);

	g4Reco->registerSubsystem(cyl);

	// Outer skin:
	cyl = new PHG4CylinderSubsystem("DIRC_SMALL_CST_Outer_Skin", 11);
	cyl->set_double_param("radius", G4DIRC_SMALL::outer_skin_radius - 0.076);
	cyl->set_double_param("length", G4DIRC_SMALL::length);
	cyl->set_string_param("material", "G4_Al");
	cyl->set_double_param("thickness", 0.076);
	cyl->set_double_param("place_x", 0.);
	cyl->set_double_param("place_y", 0.);
	cyl->set_double_param("place_z", G4DIRC_SMALL::z_shift);
	cyl->SetActive(0);
	cyl->SuperDetector("DIRC_SMALL");
	cyl->OverlapCheck(OverlapCheck);

	g4Reco->registerSubsystem(cyl);

	// Done
	return G4DIRC_SMALL::outer_skin_radius;
}
#endif
