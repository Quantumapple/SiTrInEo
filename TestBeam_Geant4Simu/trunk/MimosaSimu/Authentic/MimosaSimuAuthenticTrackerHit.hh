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
// $Id: MimosaSimuTrackerHit.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file MimosaSimuTrackerHit.hh
/// \brief Definition of the MimosaSimuTrackerHit class

#ifndef MimosaSimuAuthenticTrackerHit_h
#define MimosaSimuAuthenticTrackerHit_h 1

#include <iostream>
#include "Riostream.h"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "tls.hh"
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include "MimosaSimuTrackerHit.hh"

/// Authentic Tracker hit class

class MimosaSimuAuthenticTrackerHit : public TObject {

  // Root class derived from Tobject to store hit information

  public:
    MimosaSimuAuthenticTrackerHit();
    MimosaSimuAuthenticTrackerHit(const MimosaSimuTrackerHit&);
    MimosaSimuAuthenticTrackerHit(const MimosaSimuAuthenticTrackerHit &authenticHit);
    virtual ~MimosaSimuAuthenticTrackerHit();

  public:
    G4int           fAuthenticTrackID;
    G4int           fAuthenticTrackPDGCode;
    TVector3        fAuthenticTrackVertexPos;
    G4int           fAuthenticMimosaSensorID;
    G4double        fAuthenticEdep;
    TVector3        fAuthenticPos;
    TVector3        fAuthenticStepLength;
    TVector3        fAuthenticFirstStepPos;
    TVector3        fAuthenticLastStepPos;
    TLorentzVector  fAuthenticFirstStep4Vector;
    TLorentzVector  fAuthenticLastStep4Vector;
    
    ClassDef(MimosaSimuAuthenticTrackerHit,1)  //a MimosaSimuTrackerHit tranformed into a class derived from TObject
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
#endif
