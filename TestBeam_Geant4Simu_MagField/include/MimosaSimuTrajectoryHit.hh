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
// $Id: MimosaSimuTrajectoryHit.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file MimosaSimuTrajectoryHit.hh
/// \brief Definition of the MimosaSimuTrajectoryHit class

#ifndef MimosaSimuTrajectoryHit_h
#define MimosaSimuTrajectoryHit_h 1

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

/// Tracker hit class

class MimosaSimuTrajectoryHit : public G4VHit
{
  public:
    MimosaSimuTrajectoryHit();
    MimosaSimuTrajectoryHit(const MimosaSimuTrajectoryHit&);
    virtual ~MimosaSimuTrajectoryHit();

    // operators
    const MimosaSimuTrajectoryHit& operator=(const MimosaSimuTrajectoryHit&);
    G4int operator==(const MimosaSimuTrajectoryHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID         (G4int trackID)       { fTrackID        = trackID; };
    void SetTrackPDFCode    (G4int trackPDG)      { fTrackPDGCode   = trackPDG; };
    void SetTrackVertexPos  (G4ThreeVector xyz)   { fTrackVertexPos = xyz; };
    void SetTrackVertexTime (G4double time)       { fTrackVertexTime = time; };
    
    void SetEdep            (G4double de)         { fEdep = de; };
    void SetFirstStepPos    (G4ThreeVector xyz)   { fFirstStepPos = xyz; };
    void SetFirstStepTime   (G4double time)       { fFirstStepTime = time; };
    void SetLastStepPos     (G4ThreeVector xyz)   { fLastStepPos  = xyz; };
    void SetLastStepTime    (G4double time)       { fLastStepTime = time; };
    void SetFirstStep4Vector(G4LorentzVector exyz) { fFirstStep4Vector = exyz; };
    void SetLastStep4Vector (G4LorentzVector exyz) { fLastStep4Vector  = exyz; };

    // Get methods
    G4int GetTrackID() const                    { return fTrackID; };
    G4int GetTrackPDGCode() const               { return fTrackPDGCode; };
    G4ThreeVector GetTrackVertexPos() const     { return fTrackVertexPos; };
    G4double      GetTrackVertexTime() const    { return fTrackVertexTime; };
    
    G4double GetEdep() const                    { return fEdep; };
    G4ThreeVector GetFirstStepPos() const       { return fFirstStepPos; };
    G4double      GetFirstStepTime() const      { return fFirstStepTime; };
    G4ThreeVector GetLastStepPos() const        { return fLastStepPos; };
    G4double      GetLastStepTime() const       { return fLastStepTime; };
    G4LorentzVector GetFirstStep4Vector() const { return fFirstStep4Vector; };
    G4LorentzVector GetLastStep4Vector() const  { return fLastStep4Vector; };

  private:

      G4int           fTrackID;
      G4int           fTrackPDGCode;
      G4ThreeVector   fTrackVertexPos;
      G4double        fTrackVertexTime;
      G4double        fEdep;
      G4ThreeVector   fFirstStepPos;
      G4double        fFirstStepTime;
      G4ThreeVector   fLastStepPos;
      G4double        fLastStepTime;
      G4LorentzVector fFirstStep4Vector;
      G4LorentzVector fLastStep4Vector;
      
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MimosaSimuTrajectoryHit> MimosaSimuTrajectoryHitsCollection;

extern G4ThreadLocal G4Allocator<MimosaSimuTrajectoryHit>* MimosaSimuTrajectoryHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MimosaSimuTrajectoryHit::operator new(size_t)
{
  if(!MimosaSimuTrajectoryHitAllocator)
      MimosaSimuTrajectoryHitAllocator = new G4Allocator<MimosaSimuTrajectoryHit>;
  return (void *) MimosaSimuTrajectoryHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void MimosaSimuTrajectoryHit::operator delete(void *hit)
{
  MimosaSimuTrajectoryHitAllocator->FreeSingle((MimosaSimuTrajectoryHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
#endif
