//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: MimosaSimuDetectorConstruction.hh 73722 2013-09-09 10:23:05Z gcosmo $
//
/// \file MimosaSimuDetectorConstruction.hh
/// \brief Definition of the MimosaSimuDetectorConstruction class

#ifndef MimosaSimuDetectorConstruction_h
#define MimosaSimuDetectorConstruction_h 1

#include "globals.hh"
#include "G4GDMLParser.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"
#include "MimosaSimuSetup.hh"
#include "G4RotationMatrix.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class MimosaSimuDetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

const int      MyNumberOfSensors = 50;

class MimosaSimuDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MimosaSimuDetectorConstruction(MimosaSimuSetup* TheSetup);
    virtual ~MimosaSimuDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Set methods
    void SetMimosaSensorMaterial(G4int, G4String );
    void SetMaxStep(G4double );
    void SetCheckOverlaps(G4bool );
    void CheckMaterial(TString Material);
    
    void SetRotationMatrix(G4double RotationX,
                           G4double RotationY,
			   G4double RotationZ,
			   G4RotationMatrix* &rot);
    
  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
  
    MimosaSimuSetup* fSetup;
    
    G4double fWorldLengthX;
    G4double fWorldLengthY;
    G4double fWorldLengthZ;
    
    G4int fNumberOfSensors;
    
    //GDML parser
    G4GDMLParser fParser;
    //Name of output file with GDML geometry
    G4String fWriteFile;
    
    // data members
    G4int fNbOfChambers;

    G4LogicalVolume**    fLogicMimosaSensor;           // pointer to the logical  Mimosa Sensor (epitaxial layer => sensitive)
    G4VPhysicalVolume**  fPhysicMimosaSensor;          // pointer to the physical Mimosa Sensor (epitaxial layer => sensitive)
    
    G4LogicalVolume**    fLogicMimosaSensorMetal;      // pointer to the logical  Mimosa metal layer => non-sensitive
    G4VPhysicalVolume**  fPhysicMimosaSensorMetal;     // pointer to the physical Mimosa metal layer => non-sensitive
    
    G4LogicalVolume**    fLogicMimosaSensorSubstrate;  // pointer to the logical  Mimosa substrate layer => non-sensitive
    G4VPhysicalVolume**  fPhysicMimosaSensorSubstrate; // pointer to the physical Mimosa substrate layer => non-sensitive
    
    //Only for sources
    G4LogicalVolume*     fLogicSource;                 // pointer to source logical  volume
    G4VPhysicalVolume*   fPhysicSource;                // pointer to source physical volume

    G4Material*        fTargetMaterial;  // pointer to the target  material
    G4Material*        fChamberMaterial; // pointer to the chamber material
    G4Material**       fSensorMaterial;  // pointer to sensor material
    G4Material*        air;              // pointer to air material
    G4Material*        vacuum;           // pointer to vacuum material
    G4Material*        silicon;          // pointer to silicon material
    G4Material*        carbon;           // pointer to carbon material
    G4Material*        tape;             // pointer to tape material
    G4Material*        mylar;            // pointer to mylar material
    G4Material*        aluminium;        // pointer to aluminium material
    G4Material*        copper;           // pointer to copper material
    G4Material*        beryllium;        // pointer to beryllium material
    G4Material*        kapton;           // pointer to kapton material
    G4Material*        glue;             // pointer to glue material
    G4Material*        SiC_foam;         // pointer to Si-C foam
    

    G4UserLimits* fStepLimit;            // pointer to user step limits

    MimosaSimuDetectorMessenger* fMessenger;   // messenger

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                         // magnetic field messenger
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 
    
    //Material map
    std::map<TString,G4Material*> _material_map;
    std::vector<TString> _material_name_list;
    
    std::vector<G4String>  _ListOfSensitiveDetectors_Names;
    std::vector<G4int>     _ListOfSensitiveDetectors_ID;
    
    std::vector<G4String>          _ListOfAllLogicalVolume_Names;
    std::vector<G4LogicalVolume*>  _ListOfAllLogicalVolumes;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
