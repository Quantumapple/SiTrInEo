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
// $Id: MimosaSimuTrackerHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file MimosaSimuTrackerHit.cc
/// \brief Implementation of the MimosaSimuTrackerHit class

#include "MimosaSimuAuthenticTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

#include <iomanip>

ClassImp(MimosaSimuAuthenticTrackerHit);

MimosaSimuAuthenticTrackerHit::MimosaSimuAuthenticTrackerHit() : TObject()
{
  
  fAuthenticTrackID          = -1;
  fAuthenticTrackPDGCode     = -999;
  fAuthenticTrackVertexPos   = TVector3(0.0,0.0,0.0);
  fAuthenticMimosaSensorID   = -1;
  fAuthenticEdep             = -999.0;
  fAuthenticPos              = TVector3(0.0,0.0,0.0);
  fAuthenticStepLength       = TVector3(0.0,0.0,0.0);
  fAuthenticFirstStepPos     = TVector3(0.0,0.0,0.0);
  fAuthenticLastStepPos      = TVector3(0.0,0.0,0.0);
  fAuthenticFirstStep4Vector = TLorentzVector(0.0,0.0,0.0,0.0);
  fAuthenticLastStep4Vector  = TLorentzVector(0.0,0.0,0.0,0.0);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuAuthenticTrackerHit::MimosaSimuAuthenticTrackerHit(const MimosaSimuTrackerHit& aHit) : TObject()
{

  G4ThreeVector   tmp;
  G4LorentzVector tmp_lorentz;
  
  fAuthenticTrackID        = aHit.GetTrackID();
  fAuthenticTrackPDGCode   = aHit.GetTrackPDGCode();
  tmp = aHit.GetTrackVertexPos();
  fAuthenticTrackVertexPos = TVector3(tmp.x()/mm,tmp.y()/mm,tmp.z()/mm);
  fAuthenticMimosaSensorID = aHit.GetMimosaSensorID();
  fAuthenticEdep           = aHit.GetEdep()/GeV;
  tmp = aHit.GetPos();
  fAuthenticPos            = TVector3(tmp.x()/mm,tmp.y()/mm,tmp.z()/mm);
  tmp = aHit.GetStepLength();
  fAuthenticStepLength     = TVector3(tmp.x()/mm,tmp.y()/mm,tmp.z()/mm);
  tmp = aHit.GetFirstStepPos();
  fAuthenticFirstStepPos   = TVector3(tmp.x()/mm,tmp.y()/mm,tmp.z()/mm);
  tmp = aHit.GetLastStepPos();
  fAuthenticLastStepPos   = TVector3(tmp.x()/mm,tmp.y()/mm,tmp.z()/mm);
  tmp_lorentz = aHit.GetFirstStep4Vector();
  fAuthenticFirstStep4Vector = TLorentzVector(tmp_lorentz.px()/GeV,tmp_lorentz.py()/GeV,tmp_lorentz.pz()/GeV,tmp_lorentz.e()/GeV);
  tmp_lorentz = aHit.GetLastStep4Vector();
  fAuthenticLastStep4Vector = TLorentzVector(tmp_lorentz.px()/GeV,tmp_lorentz.py()/GeV,tmp_lorentz.pz()/GeV,tmp_lorentz.e()/GeV);

}

MimosaSimuAuthenticTrackerHit::MimosaSimuAuthenticTrackerHit(const MimosaSimuAuthenticTrackerHit &authenticHit)
{
  
  fAuthenticTrackID          = authenticHit.fAuthenticTrackID;
  fAuthenticTrackPDGCode     = authenticHit.fAuthenticTrackPDGCode;
  fAuthenticTrackVertexPos   = authenticHit.fAuthenticTrackVertexPos;
  fAuthenticMimosaSensorID   = authenticHit.fAuthenticMimosaSensorID;
  fAuthenticEdep             = authenticHit.fAuthenticEdep;
  fAuthenticPos              = authenticHit.fAuthenticPos;
  fAuthenticStepLength       = authenticHit.fAuthenticStepLength;
  fAuthenticFirstStepPos     = authenticHit.fAuthenticFirstStepPos;
  fAuthenticLastStepPos      = authenticHit.fAuthenticLastStepPos;
  fAuthenticFirstStep4Vector = authenticHit.fAuthenticFirstStep4Vector;
  fAuthenticLastStep4Vector  = authenticHit.fAuthenticLastStep4Vector;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuAuthenticTrackerHit::~MimosaSimuAuthenticTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

