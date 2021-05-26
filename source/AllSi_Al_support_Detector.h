// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ALLSI_AL_SUPPORT_DETECTOR_H
#define ALLSI_AL_SUPPORT_DETECTOR_H

#include <g4main/PHG4Detector.h>
#include <Geant4/G4Material.hh>

#include <set>
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class AllSi_Al_support_Detector : public PHG4Detector
{
	public:
		//! constructor
		AllSi_Al_support_Detector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

		//! destructor
		virtual ~AllSi_Al_support_Detector() {}

		//! construct
		void ConstructMe(G4LogicalVolume *world) override;

		void Print(const std::string &what = "ALL") const override;

		//!@name volume accessors
		//@{
		int IsInDetector(G4VPhysicalVolume *) const;
		//@}

		void SuperDetector(const std::string &name) { m_SuperDetector = name; }
		const std::string SuperDetector() const { return m_SuperDetector; }
		G4Material * element_material( std::string identifier );
		void addDetectorSection( G4LogicalVolume *logicWorld , std::string name , double sign );

	private:
		PHParameters *m_Params;

		// active volumes
		std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;

		std::string m_SuperDetector;
};

#endif // ALLSI_AL_SUPPORT_DETECTOR_H
