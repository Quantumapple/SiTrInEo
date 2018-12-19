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
// $Id: MimosaSimuTrackerSD.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file MimosaSimuTrackerSD.hh
/// \brief Definition of the MimosaSimuTrackerSD class

#ifndef MimosaSimuTrackerSD_h
#define MimosaSimuTrackerSD_h 1

#include "G4VSensitiveDetector.hh"

#include "MimosaSimuTrackerHit.hh"

#include <vector>

//class MimosaSimuDetectorConstruction;
class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// MimosaSimuTracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class MimosaSimuTrackerSD : public G4VSensitiveDetector
{
  public:
    MimosaSimuTrackerSD(const G4String& name, 
                        const G4String& hitsCollectionName,
			const G4int TheSensorID,
			const bool  verbose);
    virtual ~MimosaSimuTrackerSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);
    
    void  SetSensorID(G4int mySensorID){ sensorID = mySensorID; };
    G4int GetSensorID()                { return sensorID; };

  private:
    MimosaSimuTrackerHitsCollection* fHitsCollection;
    G4int sensorID;
    bool  verbosity;
   
    G4int           fTrackID;              // particle track ID
    G4int           fTrackPDGCode;         // Particle pdg ID producing hit
    G4int           fMimosaSensorID;       // Sendor ID
    G4double        fEdep;                 // Deposited energy by ionization
    G4ThreeVector   fPos;                  // Average (X,Y,Z) position between 1st and last step in lab frame
    G4ThreeVector   fStepLength;           // Step length
    G4ThreeVector   fFirstStepPos;         // entrance (X,Y,Z) position in lab frame
    G4ThreeVector   fFirstStepPosLocal;    // entrance (U,V,W) position in sensor local frame
    G4ThreeVector   fFirstStepLocalAngles; // entrance (theta,phi,nothing) angles in sensor local frame
    G4double        fFirstStepTime;        // time at entrance
    G4ThreeVector   fLastStepPos;          // exit (X,Y,Z) position in lab frame
    G4ThreeVector   fLastStepPosLocal;     // exit (U,V,W) position in sensor local frame
    G4double        fLastStepTime;         // time at exit
    G4LorentzVector fFirstStep4Vector;     // 4-momentum at entrance
    G4LorentzVector fLastStep4Vector;      // 4-momentum at exit
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
