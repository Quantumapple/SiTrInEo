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
// $Id: MimosaSimuEventAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file MimosaSimuEventAction.cc
/// \brief Implementation of the MimosaSimuEventAction class

#include "MimosaSimuEventAction.hh"
#include "MimosaSimuTrackerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "g4root.hh"
#include "MimosaSimuTrackerHit.hh"
#include "MimosaSimuTrajectoryHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuEventAction::MimosaSimuEventAction(MimosaSimuHistoManager* TheHisto,
					     MimosaSimuSetup*        TheSetup) : G4UserEventAction(),
                                                                                 fHisto(TheHisto),
                                                                                 fSetup(TheSetup),
                                                                                 fPrintFreq(1)
{
  
  verbosity = fSetup->GetAnalysisPar().SavePlots;
  
  fWatch.Start();
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuEventAction::~MimosaSimuEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuEventAction::BeginOfEventAction(const G4Event* event)
{

  if(verbosity) {
    G4cout << G4endl;
    G4cout << "Beging of Event " << event->GetEventID()+1 << G4endl;
  }
  
  //if(!((event->GetEventID()+1)%fPrintFreq)) G4cout << event->GetEventID()+1 << " simulated events!!!" << G4endl;
  if(!((event->GetEventID()+1)%fPrintFreq)) {
    cout << event->GetEventID()+1 << " simulated events!!!  ";
    fWatch.Print();
    fWatch.Continue();
  }
  if(fPrintFreq < 100000) {
    if((event->GetEventID()+1) == 10*fPrintFreq) fPrintFreq *= 10;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuEventAction::EndOfEventAction(const G4Event* event)
{
  
  // get number of stored trajectories
  //G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  //G4int n_trajectories = 0;
  //if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  G4int nPlanes = fSetup->TrackerParameter.Planes;
  G4int counter = 0;

  for(G4int iPlane=0;iPlane<nPlanes;iPlane++) {
    TString MyCollectionName = TString("MimosaSensorHitsCollection") + long(iPlane+1);
    for(G4int icoll=0;icoll<event->GetHCofThisEvent()->GetNumberOfCollections();icoll++) {
      TString TheCollName(event->GetHCofThisEvent()->GetHC(icoll)->GetName().data());
      if(MyCollectionName == TheCollName) {
	MimosaSimuTrackerHitsCollection* hc = (MimosaSimuTrackerHitsCollection*)event->GetHCofThisEvent()->GetHC(icoll);
        if(verbosity) G4cout << "    "  << hc->GetSize() << " hits stored in collection " << hc->GetName().data() << " for this event." << G4endl;
	for(G4int ihit=0; ihit<hc->entries();ihit++) {
          fHisto->FillMimosaArray((*hc)[ihit],counter); 
          counter++;
        }
	break;
      }
    }
  }
  
  if(verbosity) G4cout << "   total number of Mimosa Hits = " << counter << G4endl;
  
  counter = 0;
  for(G4int icoll=0;icoll<event->GetHCofThisEvent()->GetNumberOfCollections();icoll++) {
    TString TheCollName(event->GetHCofThisEvent()->GetHC(icoll)->GetName().data());
    
    bool IsSensitive = false;
    for(G4int iPlane=0;iPlane<nPlanes;iPlane++) {
      TString MyCollectionName = TString("MimosaSensorHitsCollection") + long(iPlane+1);
      if(MyCollectionName == TheCollName) {
	IsSensitive = true;
	break;
      }
    }
    
    if(IsSensitive) continue;
    
    MimosaSimuTrajectoryHitsCollection* hc = (MimosaSimuTrajectoryHitsCollection*)event->GetHCofThisEvent()->GetHC(icoll);
    if(verbosity) G4cout << "    "  << hc->GetSize() << " hits stored in collection " << hc->GetName().data() << " for this event." << G4endl;
    for(G4int ihit=0; ihit<hc->entries();ihit++) {
      fHisto->FillNonSensitiveArray((*hc)[ihit],counter); 
      counter++;
    }
    
  }
  if(verbosity) G4cout << "   total number of Non-sensitive Hits = " << counter << G4endl;
  
  fHisto->FillTreeBlocks();
  
  fHisto->FillNtuple();
  
  if(verbosity) {
    G4cout << "End of Event " << event->GetEventID()+1 << G4endl;
    G4cout << G4endl;
  }

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
