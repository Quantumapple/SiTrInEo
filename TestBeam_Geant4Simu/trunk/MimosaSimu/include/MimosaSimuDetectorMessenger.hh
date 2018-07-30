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
// $Id: MimosaSimuDetectorMessenger.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file MimosaSimuDetectorMessenger.hh
/// \brief Definition of the MimosaSimuDetectorMessenger class

#ifndef MimosaSimuDetectorMessenger_h
#define MimosaSimuDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MimosaSimuDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for MimosaSimuDetectorConstruction.
///
/// It implements commands:
/// - /MimosaSimu/det/setSensorMaterial name
/// - /MimosaSimu/det/stepMax value unit

class MimosaSimuDetectorMessenger: public G4UImessenger
{
  public:
    MimosaSimuDetectorMessenger(MimosaSimuDetectorConstruction* );
    virtual ~MimosaSimuDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MimosaSimuDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           fMimosaSimuDirectory;
    G4UIdirectory*           fDetDirectory;

    G4UIcmdWithAString*      fSensorMatCmd;

    G4UIcmdWithADoubleAndUnit* fStepMaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
