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
// $Id: MimosaSimuTrajectorySD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file MimosaSimuTrajectorySD.cc
/// \brief Implementation of the MimosaSimuTrajectorySD class

#include "MimosaSimuTrajectorySD.hh"
//#include "MimosaSimuDetectorConstruction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrajectorySD::MimosaSimuTrajectorySD(const G4String& name,
                                               const G4String& hitsCollectionName,
					       const bool  verbose) : 
                                               G4VSensitiveDetector(name),
                                               fTrajectoryHitsCollection(NULL),
                                               verbosity(verbose)
                                         
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrajectorySD::~MimosaSimuTrajectorySD() 
{
  
  //RootIO::GetInstance()->Close();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuTrajectorySD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fTrajectoryHitsCollection = new MimosaSimuTrajectoryHitsCollection(SensitiveDetectorName, collectionName[0]);
  
  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fTrajectoryHitsCollection );
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MimosaSimuTrajectorySD::ProcessHits(G4Step* aStep, 
					   G4TouchableHistory*)
{

  G4Track* aTrack                           = aStep->GetTrack();
  G4int trackID                             = aTrack->GetTrackID();
  G4int trackParentID                       = aTrack->GetParentID();
  G4ThreeVector trackMomentum               = aTrack->GetMomentum();
  G4ThreeVector trackVertexPosition         = aTrack->GetVertexPosition();
  G4double      trackVertexTime             = aTrack->GetGlobalTime();
  const G4DynamicParticle* aDynamicParticle = aTrack->GetDynamicParticle();
  G4LorentzVector a4Vector_EnergyMomentum   = aDynamicParticle->Get4Momentum();
  G4int pdgCode                             = aDynamicParticle->GetPDGcode();
  G4double Edep                             = aStep->GetTotalEnergyDeposit();
  
  G4StepPoint* preStepPoint                 = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint                = aStep->GetPostStepPoint();
  G4ThreeVector preStepPointPosition        = preStepPoint->GetPosition();
  //G4double      preStepPointTime            = preStepPoint->GetGlobalTime();
  G4double      preStepPointTime            = aTrack->GetGlobalTime() - aStep->GetDeltaTime();
  G4ThreeVector postStepPointPosition       = postStepPoint->GetPosition();
  //G4double      postStepPointTime           = postStepPoint->GetGlobalTime();
  G4double      postStepPointTime           = aTrack->GetGlobalTime();

  G4LorentzVector preStep4Vector            = G4LorentzVector(preStepPoint->GetMomentum().x(),
                                                              preStepPoint->GetMomentum().y(),
							      preStepPoint->GetMomentum().z(),
							      preStepPoint->GetTotalEnergy());
  G4LorentzVector postStep4Vector           = G4LorentzVector(postStepPoint->GetMomentum().x(),
                                                              postStepPoint->GetMomentum().y(),
							      postStepPoint->GetMomentum().z(),
							      postStepPoint->GetTotalEnergy());
  
#if 0
  G4cout << endl;
  G4cout << "trkID = " << trackID << "  " << G4endl
         << "PDG   = " << pdgCode << "  " << G4endl
	 << "VtxPos = (" << trackVertexPosition[0] << "," << trackVertexPosition[1] << "," << trackVertexPosition[2] << ") mm  " << G4endl
	 << "Edep   = "  << Edep  << "  " << G4endl
	 << "PreStepPos[(x,y,z);t]     = [(" << preStepPointPosition[0]  << "," << preStepPointPosition[1]  << "," << preStepPointPosition[2]  << ")mm; " << preStepPointTime  << " ns]  " << G4endl
	 << "PostStepPos[(x,y,z);t]    = [(" << postStepPointPosition[0] << "," << postStepPointPosition[1] << "," << postStepPointPosition[2] << ")mm; " << postStepPointTime << " ns]  " << G4endl
	 << "PreStep4Vect(Px,Py,Pz,E)  = (" << preStep4Vector.px() << "," << preStep4Vector.py() << "," << preStep4Vector.pz() << "," << preStep4Vector.e() << ") MeV " << G4endl
	 << "PostStep4Vect(Px,Py,Pz,E) = (" << postStep4Vector.px() << "," << postStep4Vector.py() << "," << postStep4Vector.pz() << "," << postStep4Vector.e() << ") MeV " << G4endl
         << G4endl;
  G4cout << endl;
#endif
  
  MimosaSimuTrajectoryHit* myHit = new MimosaSimuTrajectoryHit();
  myHit->SetTrackID(trackID);
  myHit->SetTrackPDFCode(pdgCode);
  myHit->SetTrackVertexPos(trackVertexPosition);
  myHit->SetTrackVertexTime(trackVertexTime);
  myHit->SetEdep(Edep);
  myHit->SetFirstStepPos(preStepPointPosition);
  myHit->SetFirstStepTime(preStepPointTime);
  myHit->SetLastStepPos(postStepPointPosition);
  myHit->SetLastStepTime(postStepPointTime);
  myHit->SetFirstStep4Vector(preStep4Vector);
  myHit->SetLastStep4Vector(postStep4Vector);
     
  fTrajectoryHitsCollection->insert(myHit);

  return true;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuTrajectorySD::EndOfEvent(G4HCofThisEvent*)
{
  
  G4int nofHits = fTrajectoryHitsCollection->entries();
  if(verbosity) {
    for(G4int i=0; i<nofHits; i++ ) (*fTrajectoryHitsCollection)[i]->Print();
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
