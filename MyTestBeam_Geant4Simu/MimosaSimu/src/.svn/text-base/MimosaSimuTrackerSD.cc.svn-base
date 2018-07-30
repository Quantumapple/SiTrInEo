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
// $Id: MimosaSimuTrackerSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file MimosaSimuTrackerSD.cc
/// \brief Implementation of the MimosaSimuTrackerSD class

#include "MimosaSimuTrackerSD.hh"
//#include "MimosaSimuDetectorConstruction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "g4root.hh"
#include "G4TouchableHandle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrackerSD::MimosaSimuTrackerSD(const G4String& name,
                                         const G4String& hitsCollectionName,
					 const G4int TheSensorID,
					 const bool  verbose) : 
                                         G4VSensitiveDetector(name),
                                         fHitsCollection(NULL),
                                         sensorID(TheSensorID),
                                         verbosity(verbose)
                                         
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuTrackerSD::~MimosaSimuTrackerSD() 
{
  
  //RootIO::GetInstance()->Close();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection = new MimosaSimuTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  
  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MimosaSimuTrackerSD::ProcessHits(G4Step* aStep, 
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

  G4StepPoint* preStepPoint  = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();

  G4ThreeVector preStepPointPosition    = preStepPoint->GetPosition();
  G4double      preStepPointGlobalTime  = trackVertexTime - aStep->GetDeltaTime();
  G4ThreeVector postStepPointPosition   = postStepPoint->GetPosition();
  G4double      postStepPointGlobalTime = trackVertexTime;

  if(preStepPoint->GetStepStatus() == fGeomBoundary) {
    //G4cout << "First Step in Hit" << G4endl;
    fEdep = 0.0;
    fEdep    += aStep->GetTotalEnergyDeposit();
    fFirstStepPos     = preStepPointPosition;
    fFirstStepTime    = preStepPointGlobalTime;
    fFirstStep4Vector = a4Vector_EnergyMomentum;
    
    // We want the track angles of incidence on the local volume and the local position
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4ThreeVector worldDirection   = preStepPoint->GetMomentumDirection();
    G4ThreeVector localDirection   = theTouchable->GetHistory()->GetTopTransform().TransformAxis(worldDirection);
   
    // Calculating the local particle's incidence angles
    G4double theta_local,phi_local;
    theta_local = acos(localDirection.z());
    
    G4double epsilon = 1.0e-10;
    if(TMath::Abs(localDirection.x()) < epsilon) {
      if(TMath::Abs(localDirection.y()) < epsilon) phi_local =  0.0;
      else if(localDirection.y() > 0.0)            phi_local =  TMath::Pi()/2.0;
      else if(localDirection.y() < 0.0)            phi_local = -TMath::Pi()/2.0;
    }
    else {
      G4double phi_temp = atan(TMath::Abs(localDirection.y()/localDirection.x()));
      if(localDirection.x() > 0) {
	if(TMath::Abs(localDirection.y()) < epsilon) phi_local =  0.0;
	else if(localDirection.y() > 0.0)            phi_local =  phi_temp;
	else if(localDirection.y() < 0.0)            phi_local = -phi_temp;
      }
      else if(localDirection.x() < 0) {
	if(TMath::Abs(localDirection.y()) < epsilon) phi_local =  TMath::Pi();
	else if(localDirection.y() > 0.0)            phi_local =  TMath::Pi() - phi_temp;
	else if(localDirection.y() < 0.0)            phi_local = -TMath::Pi() + phi_temp;
      }
    }
    fFirstStepLocalAngles = G4ThreeVector(theta_local,phi_local,0.0);
    //cout << "theta_local = " << theta_local/deg << ", Z-direction     = "  << localDirection.z() << endl;
    //cout << "phi_local   = " << phi_local/deg   << ", (X,Y)-direction = (" << localDirection.x() << "," << localDirection.y() << ")" << endl;

    G4ThreeVector localINPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(preStepPointPosition);
    fFirstStepPosLocal.setX(localINPosition.x());
    fFirstStepPosLocal.setY(localINPosition.y());
    fFirstStepPosLocal.setZ(localINPosition.z());
  }
  if(preStepPoint->GetStepStatus() != fGeomBoundary && postStepPoint->GetStepStatus() != fGeomBoundary) {
    //G4cout << "Intermediate Step in Hit" << G4endl;
    fEdep += aStep->GetTotalEnergyDeposit();
  }
 
  if(postStepPoint->GetStepStatus() == fGeomBoundary || aTrack->GetTrackStatus() != fAlive) {
    if(aTrack->GetTrackStatus() != fAlive && verbosity) G4cout << "The track status is different from fAlive!!!" << G4endl;
    if(preStepPoint->GetStepStatus() != fGeomBoundary) fEdep += aStep->GetTotalEnergyDeposit();
    //G4cout << "Last Step in Hit" << G4endl;
    
    fLastStepPos  = postStepPointPosition;
    fLastStepTime = postStepPointGlobalTime;
    fStepLength   = G4ThreeVector(fLastStepPos.x() - fFirstStepPos.x(),
                                  fLastStepPos.y() - fFirstStepPos.y(),
				  fLastStepPos.z() - fFirstStepPos.z());
    
    G4ThreeVector localOUTPosition = preStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(postStepPointPosition);
    fLastStepPosLocal.setX(localOUTPosition.x());
    fLastStepPosLocal.setY(localOUTPosition.y());
    fLastStepPosLocal.setZ(localOUTPosition.z());
    
    fLastStep4Vector = a4Vector_EnergyMomentum;
    
    G4ThreeVector Position(0.5*(fFirstStepPos.x() + fLastStepPos.x()),
                           0.5*(fFirstStepPos.y() + fLastStepPos.y()),
			   0.5*(fFirstStepPos.z() + fLastStepPos.z()));
    
    if(verbosity) {
      G4cout << "Sensor ID   = "  << GetSensorID() << G4endl;
      G4cout << "PreStepPos  = (" << fFirstStepPos[0]/um << "," << fFirstStepPos[1]/um << "," << fFirstStepPos[2]/um << ") um" << G4endl;
      G4cout << "PostStepPos = (" << fLastStepPos[0]/um  << "," << fLastStepPos[1]/um  << "," << fLastStepPos[2]/um  << ") um" << G4endl;
    }
    
#if 0
    cout << "ProcessHits:: trackID = " << trackID << ", "
         << "Momentum = " << sqrt(pow(fFirstStep4Vector.px(),2) + pow(fFirstStep4Vector.py(),2) + pow(fFirstStep4Vector.pz(),2))/MeV << " MeV/c,  "
         << "pdgID = " << pdgCode << ",  "
	 << "sensorID = " << GetSensorID() << ",  "
         << endl;
#endif

    MimosaSimuTrackerHit* myHit = new MimosaSimuTrackerHit();
    myHit->SetTrackID(trackID);
    myHit->SetTrackPDFCode(pdgCode);
    myHit->SetTrackVertexPos(trackVertexPosition);
    myHit->SetTrackVertexTime(trackVertexTime);
    myHit->SetMimosaSensorID(GetSensorID());
    myHit->SetEdep(fEdep);
    myHit->SetPos(Position);
    myHit->SetStepLength(fStepLength);
    myHit->SetFirstStepPos(fFirstStepPos);
    myHit->SetFirstStepPosLocal(fFirstStepPosLocal);
    myHit->SetFirstStepLocalAngles(fFirstStepLocalAngles);
    myHit->SetFirstStepTime(fFirstStepTime);
    myHit->SetLastStepPos(fLastStepPos);
    myHit->SetLastStepPosLocal(fLastStepPosLocal);
    myHit->SetLastStepTime(fLastStepTime);
    myHit->SetFirstStep4Vector(fFirstStep4Vector);
    myHit->SetLastStep4Vector(fLastStep4Vector);
    
    fHitsCollection->insert(myHit);
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  
  G4int nofHits = fHitsCollection->entries();
  if(verbosity) {
    for(G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
