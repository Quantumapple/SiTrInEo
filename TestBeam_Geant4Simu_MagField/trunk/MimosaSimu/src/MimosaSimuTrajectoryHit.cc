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
// $Id: MimosaSimuTrajectoryHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file MimosaSimuTrajectoryHit.cc
/// \brief Implementation of the MimosaSimuTrajectoryHit class

#include "MimosaSimuTrajectoryHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<MimosaSimuTrajectoryHit>* MimosaSimuTrajectoryHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrajectoryHit::MimosaSimuTrajectoryHit()
 : G4VHit(),
   fTrackID(-1),
   fTrackPDGCode(-1),
   fTrackVertexPos(G4ThreeVector()),
   fTrackVertexTime(-1),
   fEdep(0.),
   fFirstStepPos(G4ThreeVector()),
   fFirstStepTime(0.0),
   fLastStepPos(G4ThreeVector()),
   fLastStepTime(0.0),
   fFirstStep4Vector(G4LorentzVector(0.0,0.0,0.0,0.0)),
   fLastStep4Vector(G4LorentzVector(0.0,0.0,0.0,0.0))
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrajectoryHit::~MimosaSimuTrajectoryHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrajectoryHit::MimosaSimuTrajectoryHit(const MimosaSimuTrajectoryHit& right)
  : G4VHit()
{
  fTrackID          = right.fTrackID;
  fTrackPDGCode     = right.fTrackPDGCode;
  fTrackVertexPos   = right.fTrackVertexPos;
  fTrackVertexTime  = right.fTrackVertexTime;
  
  fEdep             = right.fEdep;
  fFirstStepPos     = right.fFirstStepPos;
  fFirstStepTime    = right.fFirstStepTime;
  fLastStepPos      = right.fLastStepPos;
  fLastStepTime     = right.fLastStepTime;
  fFirstStep4Vector = right.fFirstStep4Vector;
  fLastStep4Vector  = right.fLastStep4Vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const MimosaSimuTrajectoryHit& MimosaSimuTrajectoryHit::operator=(const MimosaSimuTrajectoryHit& right)
{
  fTrackID          = right.fTrackID;
  fTrackPDGCode     = right.fTrackPDGCode;
  fTrackVertexPos   = right.fTrackVertexPos;
  fTrackVertexTime  = right.fTrackVertexTime;
  
  fEdep             = right.fEdep;
  fFirstStepPos     = right.fFirstStepPos;
  fFirstStepTime    = right.fFirstStepTime;
  fLastStepPos      = right.fLastStepPos;
  fLastStepTime     = right.fLastStepTime;
  fFirstStep4Vector = right.fFirstStep4Vector;
  fLastStep4Vector  = right.fLastStep4Vector;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MimosaSimuTrajectoryHit::operator==(const MimosaSimuTrajectoryHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuTrajectoryHit::Draw()
{
#if 0
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
    
    G4Circle circleFirst(fFirstStepPos);
    circleFirst.SetScreenSize(4.);
    circleFirst.SetFillStyle(G4Circle::filled);
    colour = G4Colour(0.,0.,1.);
    G4VisAttributes attribsFirst(colour);
    circleFirst.SetVisAttributes(attribsFirst);
    pVVisManager->Draw(circleFirst);
    
    G4Circle circleLast(fLastStepPos);
    circleLast.SetScreenSize(4.);
    circleLast.SetFillStyle(G4Circle::filled);
    colour = G4Colour(0.,1.,0.);
    G4VisAttributes attribsLast(colour);
    circleLast.SetVisAttributes(attribsLast);
    pVVisManager->Draw(circleLast);
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuTrajectoryHit::Print()
{
  
  G4int Ndigits = 5;
  
  G4cout
     << "  trackID: " << fTrackID << G4endl
     << "  trackPDG: " << fTrackPDGCode << G4endl
     << "  trackVertexPos[(x,y,z);y]: (" << fTrackVertexPos.x()/cm << "," << fTrackVertexPos.y()/cm << "," << fTrackVertexPos.z()/cm << ") cm; " << fTrackVertexTime/ns << " ns]" << G4endl
     << "  Edep:                            "
     << std::setw(Ndigits) << fEdep/keV << " keV" << G4endl
     << "  1st- step-Position[(x,y,z);t]:   [("
     << std::setw(Ndigits) << fFirstStepPos.x()/cm << "," << fFirstStepPos.y()/cm << "," << fFirstStepPos.z()/cm << ") cm; " << fFirstStepTime/ns << " ns]" << G4endl
     << "  last-step-Position[(x,y,z);t]:   [("
     << std::setw(Ndigits) << fLastStepPos.x()/cm << "," << fLastStepPos.y()/cm << "," << fLastStepPos.z()/cm << ") cm; " << fLastStepTime/ns << " ns]" << G4endl
     << "  1st -step-4Momentum(px,py,pz,E): ("
     << std::setw(Ndigits) << fFirstStep4Vector.px()/GeV << "," << fFirstStep4Vector.py()/GeV << "," << fFirstStep4Vector.pz()/GeV << "," << fFirstStep4Vector.e()/GeV << ") GeV" << G4endl
     << "  last-step-4Momentum(px,py,pz,E): ("
     << std::setw(Ndigits) << fLastStep4Vector.px()/GeV << "," << fLastStep4Vector.py()/GeV << "," << fLastStep4Vector.pz()/GeV << "," << fLastStep4Vector.e()/GeV << ") GeV" << G4endl
     << G4endl;
     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

