#include "AllSi_vtx_serv_2lyr_Detector.h"
#include <phparameter/PHParameters.h>
#include <g4main/PHG4Detector.h>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;
using namespace std;
//____________________________________________________________________________..
AllSi_vtx_serv_2lyr_Detector::AllSi_vtx_serv_2lyr_Detector(PHG4Subsystem *subsys,
		PHCompositeNode *Node,
		PHParameters *parameters,
		const std::string &dnam)
	: PHG4Detector(subsys, Node, dnam)
	  , m_Params(parameters)
{
}
//_______________________________________________________________
int AllSi_vtx_serv_2lyr_Detector::IsInDetector(G4VPhysicalVolume *volume) const
{
	set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if (iter != m_PhysicalVolumesSet.end())
	{
		return 1;
	}
	return 0;
}
//_______________________________________________________________
void AllSi_vtx_serv_2lyr_Detector::ConstructMe(G4LogicalVolume *logicWorld)
{
	//begin implement your own here://
	// Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
	addDetectorSection( logicWorld , "Al_pos_z" ,  1. );
	addDetectorSection( logicWorld , "Al_neg_z" , -1. );
	//end implement your own here://
	return;
}
// ======================================================================================================
void AllSi_vtx_serv_2lyr_Detector::Print(const std::string &what) const
{
	cout << "AllSi_vtx_serv_2lyr_ Detector:" << endl;
	if (what == "ALL" || what == "VOLUME")
	{
		cout << "Version 0.1" << endl;
		cout << "Parameters:" << endl;
		m_Params->Print();
	}
	return;
}
// ======================================================================================================
void  AllSi_vtx_serv_2lyr_Detector::addDetectorSection( G4LogicalVolume *logicWorld , std::string name , double sign){

	//double z_det[3] = {20.,43.23*25./18.5, 121.};
	//double rin  [3] = {18.5 * 20.0 / 25.0,43.23,43.23};
	//double rout [3] = {0};
	
	// Second vertexing layer
	double z_det_l2[2] = {15.,20.};
   	double rin_l2  [2] = {5.48,14.8};
	double rout_l2 [2] = {0};
	
	// First vertexing layer
	double z_det_l1[2] = {15.,20.};
        double rin_l1  [2] = {3.3,14.36};
        double rout_l1 [2] = {0};
	
	const int nzplanes_l2 = sizeof(z_det_l2)/sizeof(*z_det_l2);
	for(int i = 0 ; i < nzplanes_l2 ; i++){
		rout_l2 [i] = rin_l2[i]+0.44;
		z_det_l2[i] *= 10.*sign/abs(sign);
		rin_l2  [i] *= 10.;
		rout_l2 [i] *= 10.;
	}
	
	const int nzplanes_l1 = sizeof(z_det_l1)/sizeof(*z_det_l1);
        for(int i = 0 ; i < nzplanes_l1 ; i++){	
		rout_l1 [i] = rin_l1[i]+0.44;
		z_det_l1[i] *= 10.*sign/abs(sign);
                rin_l1  [i] *= 10.;
                rout_l1 [i] *= 10.;
	}
	G4Material * G4_mat_l2 = G4Material::GetMaterial("G4_Al");
	G4Material * G4_mat_l1 = G4Material::GetMaterial("G4_GRAPHITE");

	G4RotationMatrix *rotm = new G4RotationMatrix();
	rotm->rotateX(0);
	rotm->rotateY(0);
	rotm->rotateZ(0);

	G4Color G4_color_l2 = G4Color(G4Colour::Yellow());
	G4Color G4_color_l1 = G4Color(G4Colour::Cyan());

	G4VSolid *G4_polycone_l2 = new G4Polycone(name+"_l2",0,360*degree,nzplanes_l2,z_det_l2,rin_l2,rout_l2);
	G4LogicalVolume *logical_l2 = new G4LogicalVolume(G4_polycone_l2,G4_mat_l2, "AllSi_vtx_serv_2lyr_Logical");
	G4VisAttributes *vis_l2 = new G4VisAttributes(G4_color_l2);
	vis_l2->SetForceSolid(true);
	logical_l2->SetVisAttributes(vis_l2);

	G4VSolid *G4_polycone_l1 = new G4Polycone(name+"_l1",0,360*degree,nzplanes_l1,z_det_l1,rin_l1,rout_l1);
        G4LogicalVolume *logical_l1 = new G4LogicalVolume(G4_polycone_l1,G4_mat_l1, "AllSi_vtx_serv_2lyr_Logical");
        G4VisAttributes *vis_l1 = new G4VisAttributes(G4_color_l1);
        vis_l1->SetForceSolid(true);
	logical_l1->SetVisAttributes(vis_l1);

	G4VPhysicalVolume *phy_l2 = new G4PVPlacement(rotm, G4ThreeVector(0,0,0), logical_l2 , "AllSi_vtx_serv_2lyr_", logicWorld, 0, false, OverlapCheck());
	G4VPhysicalVolume *phy_l1 = new G4PVPlacement(rotm, G4ThreeVector(0,0,0), logical_l1 , "AllSi_vtx_serv_2lyr_", logicWorld, 0, false, OverlapCheck());

	// add it to the list of placed volumes so the IsInDetector method picks them up
	m_PhysicalVolumesSet.insert(phy_l2);
	m_PhysicalVolumesSet.insert(phy_l1);
}
