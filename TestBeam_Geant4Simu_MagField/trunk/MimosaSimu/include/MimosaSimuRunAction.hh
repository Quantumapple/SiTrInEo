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
// $Id: MimosaSimuRunAction.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MimosaSimuRunAction.hh
/// \brief Definition of the MimosaSimuRunAction class

#ifndef MimosaSimuRunAction_h
#define MimosaSimuRunAction_h 1

#include "g4root.hh"
#include "MimosaSimuSetup.hh"
#include "MimosaSimuHistoManager.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

/// Run action class

class MimosaSimuRunAction : public G4UserRunAction
{
  public:
    MimosaSimuRunAction(MimosaSimuSetup* TheSetup,
                        MimosaSimuHistoManager* TheHisto);
    virtual ~MimosaSimuRunAction();

    virtual void  BeginOfRunAction(const G4Run* run);
    virtual void  EndOfRunAction(const G4Run* run);
    
  private:
    MimosaSimuSetup* fSetup;
    MimosaSimuHistoManager* fHisto;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
