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
// $Id: MimosaSimuDetectorMessenger.cc 69706 2013-05-13 09:12:40Z gcosmo $
// 
/// \file MimosaSimuDetectorMessenger.cc
/// \brief Implementation of the MimosaSimuDetectorMessenger class

#include "MimosaSimuDetectorMessenger.hh"
#include "MimosaSimuDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuDetectorMessenger::MimosaSimuDetectorMessenger(MimosaSimuDetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fMimosaSimuDirectory = new G4UIdirectory("/MimosaSimu/");
  fMimosaSimuDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/MimosaSimu/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fSensorMatCmd = new G4UIcmdWithAString("/MimosaSimu/det/setSensorMaterial",this);
  fSensorMatCmd->SetGuidance("Select Material of the Sensor.");
  fSensorMatCmd->SetParameterName("choice",false);
  fSensorMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/MimosaSimu/det/stepMax",this);
  fStepMaxCmd->SetGuidance("Define a step max");
  fStepMaxCmd->SetParameterName("stepMax",false);
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuDetectorMessenger::~MimosaSimuDetectorMessenger()
{
  delete fSensorMatCmd;
  delete fStepMaxCmd;
  delete fMimosaSimuDirectory;
  delete fDetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if( command == fSensorMatCmd) {
    fDetectorConstruction->SetMimosaSensorMaterial(0,newValue);
  }

  if( command == fStepMaxCmd ) {
    fDetectorConstruction->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
