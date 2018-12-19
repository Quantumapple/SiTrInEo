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
// $Id: MimosaSimuRunAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file MimosaSimuRunAction.cc
/// \brief Implementation of the MimosaSimuRunAction class

#include "MimosaSimuRunAction.hh"

#include "TString.h"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuRunAction::MimosaSimuRunAction(MimosaSimuSetup* TheSetup,
                                         MimosaSimuHistoManager* TheHisto) : G4UserRunAction(),
                                                                             fSetup(TheSetup),
                                                                             fHisto(TheHisto)
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuRunAction::~MimosaSimuRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuRunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  //Booking the n-tuples and histos
  fHisto->book();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuRunAction::EndOfRunAction(const G4Run* )
{
  //Saving the n-tuples and histos
  fHisto->save(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
