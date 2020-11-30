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
// $Id: MimosaSimuEventAction.hh 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file MimosaSimuEventAction.hh
/// \brief Definition of the MimosaSimuEventAction class

#ifndef MimosaSimuEventAction_h
#define MimosaSimuEventAction_h 1

#include "G4UserEventAction.hh"
#include "MimosaSimuHistoManager.hh"
#include "MimosaSimuSetup.hh"

#include "TStopwatch.h"

#include "globals.hh"

/// Event action class

class MimosaSimuEventAction : public G4UserEventAction {
  public:
    MimosaSimuEventAction(MimosaSimuHistoManager *TheHisto, MimosaSimuSetup *TheSetup);
    virtual ~MimosaSimuEventAction();

    virtual void BeginOfEventAction(const G4Event *);
    virtual void EndOfEventAction(const G4Event *);

  private:
    MimosaSimuHistoManager *fHisto;
    MimosaSimuSetup *fSetup;
    bool verbosity;
    G4int fPrintFreq;
    TStopwatch fWatch;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
