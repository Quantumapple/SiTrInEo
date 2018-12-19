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

#ifndef MimosaSimuTrackerHit_h
#define MimosaSimuTrackerHit_h 1

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

class MimosaSimuTrackerHit : public G4VHit
{
  public:
    MimosaSimuTrackerHit();
    MimosaSimuTrackerHit(const MimosaSimuTrackerHit&);
    virtual ~MimosaSimuTrackerHit();

    // operators
    const MimosaSimuTrackerHit& operator=(const MimosaSimuTrackerHit&);
    G4int operator==(const MimosaSimuTrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID               (G4int trackID)        { fTrackID              = trackID; };
    void SetTrackPDFCode          (G4int trackPDG)       { fTrackPDGCode         = trackPDG; };
    void SetTrackVertexPos        (G4ThreeVector xyz)    { fTrackVertexPos       = xyz; };
    void SetTrackVertexTime       (G4double time)        { fTrackVertexTime      = time; };
    
    void SetMimosaSensorID        (G4int sensorID)       { fMimosaSensorID       = sensorID; };
    void SetEdep                  (G4double de)          { fEdep                 = de; };
    void SetPos                   (G4ThreeVector xyz)    { fPos                  = xyz; };
    void SetStepLength            (G4ThreeVector xyz)    { fStepLength           = xyz; };
    void SetFirstStepPos          (G4ThreeVector xyz)    { fFirstStepPos         = xyz; };
    void SetFirstStepPosLocal     (G4ThreeVector xyz)    { fFirstStepPosLocal    = xyz; };
    void SetFirstStepLocalAngles  (G4ThreeVector xyz)    { fFirstStepLocalAngles = xyz; };
    void SetFirstStepTime         (G4double time)        { fFirstStepTime        = time; };
    void SetLastStepPos           (G4ThreeVector xyz)    { fLastStepPos          = xyz; };
    void SetLastStepPosLocal      (G4ThreeVector xyz)    { fLastStepPosLocal     = xyz; };
    void SetLastStepTime          (G4double time)        { fLastStepTime         = time; };
    void SetFirstStep4Vector      (G4LorentzVector exyz) { fFirstStep4Vector     = exyz; };
    void SetLastStep4Vector       (G4LorentzVector exyz) { fLastStep4Vector      = exyz; };

    // Get methods
    G4int GetTrackID() const                       { return fTrackID; };
    G4int GetTrackPDGCode() const                  { return fTrackPDGCode; };
    G4ThreeVector GetTrackVertexPos() const        { return fTrackVertexPos; };
    G4double      GetTrackVertexTime() const       { return fTrackVertexTime; };
    
    G4int GetMimosaSensorID() const                { return fMimosaSensorID; };
    G4double GetEdep() const                       { return fEdep; };
    G4ThreeVector GetPos() const                   { return fPos; };
    G4ThreeVector GetStepLength() const            { return fStepLength; };
    G4ThreeVector GetFirstStepPos() const          { return fFirstStepPos; };
    G4ThreeVector GetFirstStepPosLocal() const     { return fFirstStepPosLocal; };
    G4ThreeVector GetFirstStepLocalAngles() const  { return fFirstStepLocalAngles; };
    G4double      GetFirstStepTime() const         { return fFirstStepTime; };
    G4ThreeVector GetLastStepPos() const           { return fLastStepPos; };
    G4ThreeVector GetLastStepPosLocal() const      { return fLastStepPosLocal; };
    G4double      GetLastStepTime() const          { return fLastStepTime; };
    G4LorentzVector GetFirstStep4Vector() const    { return fFirstStep4Vector; };
    G4LorentzVector GetLastStep4Vector() const     { return fLastStep4Vector; };

  private:

      G4int           fTrackID;              // particle track ID
      G4int           fTrackPDGCode;         // Particle pdg ID producing hit
      G4ThreeVector   fTrackVertexPos;       // Track vertex position
      G4double        fTrackVertexTime;      // Track vertex time
      G4int           fMimosaSensorID;       // SensorID
      G4double        fEdep;                 // Deposited energy by ionization
      G4ThreeVector   fPos;                  // Average (X,Y,Z) position between 1st and last step in lab frame
      G4ThreeVector   fStepLength;           // Step length
      G4ThreeVector   fFirstStepPos;         // entrance (X,Y,Z) position in lab frame
      G4ThreeVector   fFirstStepPosLocal;    // entrance (U,V,W) position in sensor local frame
      G4ThreeVector   fFirstStepLocalAngles; // (theta,phi,nothing) in sensor local frame
      G4double        fFirstStepTime;        // time at entrance
      G4ThreeVector   fLastStepPos;          // exit (X,Y,Z) position in lab frame
      G4ThreeVector   fLastStepPosLocal;     // exit (U,V,W) position in sensor local frame
      G4double        fLastStepTime;         // time at exit
      G4LorentzVector fFirstStep4Vector;     // 4-momentum at entrance
      G4LorentzVector fLastStep4Vector;      // 4-momentum at exit
      
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MimosaSimuTrackerHit> MimosaSimuTrackerHitsCollection;

extern G4ThreadLocal G4Allocator<MimosaSimuTrackerHit>* MimosaSimuTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MimosaSimuTrackerHit::operator new(size_t)
{
  if(!MimosaSimuTrackerHitAllocator)
      MimosaSimuTrackerHitAllocator = new G4Allocator<MimosaSimuTrackerHit>;
  return (void *) MimosaSimuTrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void MimosaSimuTrackerHit::operator delete(void *hit)
{
  MimosaSimuTrackerHitAllocator->FreeSingle((MimosaSimuTrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
#endif
