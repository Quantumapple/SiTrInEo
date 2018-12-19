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
/// \file analysis/AnaEx02/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 74272 2013-/*1*/0-02 14:48:50Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MimosaSimuHistoManager.hh"
#include "TClonesArray.h"
#include "TCanvas.h"

#include <assert.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuHistoManager::MimosaSimuHistoManager(MimosaSimuSetup* TheSetup) : fRootFile(0), 
                                                                            fNtuple(0)
{
  
  fSetup = TheSetup;
  
  // ntuple
  fNtuple = 0;
  
  //Distances are in mm, Momentum/Energy in GeV, time in ns, angles in deg
  DistFactor  = mm;
  AngleFactor = deg;
  EPFactor    = GeV;
  TimeFactor  = ns;
  Edep_correction_factor = 80.0/65.0;             // possible correction factor to deposited energy
  
  EvtNumber   = 0; //Initialize event number counter
  
  //Initialization of sensitive planes (the ones with Readout > 0)
  InitChargedParticleList();
  InitSensitivePlaneList();

  //Clear up some lists
  HitsListInMimosa.clear();               // Hits in mimosa (sensitive volumes)
  HitsListInNonSensitive.clear();         // Hits in non-sensitive volumes
  ListOfParticlesInMimosa.clear();        // List of particles hitting mimosas (sensitive volumes)
  ListOfParticlesInNonSensitive.clear();  // List of particles hitting the non-sensitive volumes  
  ListOfPixelsFromNoise.clear();          // List of pixels from noise, all sensitive planes (Readout > 0)
  
  FillNonSensitiveBranch = fSetup->TrackerParameter.FillNonSensitiveBranch;  // bool used to fill and save information of hits in non-sensitive volumes
  
  rand = new TRandom(fSetup->GetAnalysisPar().MCSeed + 1809442);             // intialization of random number generator

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuHistoManager::~MimosaSimuHistoManager()
{
  
  // Deleting pointers and clearing up listss

  if(fRootFile) {
    delete fRootFile;
    fRootFile = NULL;
  }
  if(rand) {
    delete rand;
    rand = NULL;
  }
  
  HitsListInMimosa.clear();
  HitsListInNonSensitive.clear();
  
  ListOfParticlesInMimosa.clear();
  ListOfParticlesInNonSensitive.clear();
  
  for(int iplane=0;iplane<int(ListOfSensitivePlanes.size());iplane++) {
    if(ListOfSensitivePlanes[iplane].digit) {
      delete  ListOfSensitivePlanes[iplane].digit;
      ListOfSensitivePlanes[iplane].digit = NULL;
    }
  }
  ListOfSensitivePlanes.clear();
  
  ListOfPixelsFromNoise.clear();

  //cout << "MimosaSimuHistoManager::~MimosaSimuHistoManager:: All OK!!!" << endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::book()
{ 
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  TString ROOTName = fSetup->GetConfigFileName() + TString(".root");
  G4String fileName = ROOTName.Data();
  fRootFile = new TFile(fileName,"RECREATE");
  if(!fRootFile) {
    G4cout << " MimosaSimuHistoManager::book :" 
           << " problem creating the ROOT TFile "
           << G4endl;
    return;
  }

  //create ntuple
 
  fNtuple = new TTree("ntp1", "Hits");
  fNtuple->SetAutoSave(100000000);

  //==============================
  //    Mimosa hits branches
  //==============================
  //Definition of Particles Branch
  fNtuple->Branch("ParticleNb",         &ParticleNb,          "ParticleNb/I");
  fNtuple->Branch("ParticlepdgID",       ParticlepdgID,       "ParticlepdgID[ParticleNb]/I");
  fNtuple->Branch("ParticleBKGType",     ParticleBKGType,     "ParticleBKGType[ParticleNb]/I");
  fNtuple->Branch("ParticleTrkVtX",      ParticleTrkVtX,      "ParticleTrkVtX[ParticleNb]/F");
  fNtuple->Branch("ParticleTrkVtY",      ParticleTrkVtY,      "ParticleTrkVtY[ParticleNb]/F");
  fNtuple->Branch("ParticleTrkVtZ",      ParticleTrkVtZ,      "ParticleTrkVtZ[ParticleNb]/F");
  fNtuple->Branch("ParticleNHits",       ParticleNHits,       "ParticleNHits[ParticleNb]/I");
  fNtuple->Branch("ParticleFirstHitIdx", ParticleFirstHitIdx, "ParticleFirstHitIdx[ParticleNb]/I");
   
  //Definition of Hits Branch
  fNtuple->Branch("HitNb",              &HitNb,               "HitNb/I");
  fNtuple->Branch("HitParticleIdx",      HitParticleIdx,      "HitParticleIdx[HitNb]/I");
  fNtuple->Branch("HitsensorID",         HitsensorID,         "HitsensorID[HitNb]/I");
  fNtuple->Branch("HitladderID",         HitladderID,         "HitladderID[HitNb]/I");
  fNtuple->Branch("HitmoduleID",         HitmoduleID,         "HitmoduleID[HitNb]/I");
  fNtuple->Branch("HitposINmmX",         HitposINmmX,         "HitposINmmX[HitNb]/F");
  fNtuple->Branch("HitposINmmY",         HitposINmmY,         "HitposINmmY[HitNb]/F");
  fNtuple->Branch("HitposINmmZ",         HitposINmmZ,         "HitposINmmZ[HitNb]/F");
  fNtuple->Branch("HitposAVRmmU",        HitposAVRmmU,        "HitposAVRmmU[HitNb]/F");
  fNtuple->Branch("HitposAVRmmV",        HitposAVRmmV,        "HitposAVRmmV[HitNb]/F");
  fNtuple->Branch("HitposAVRmmW",        HitposAVRmmW,        "HitposAVRmmW[HitNb]/F");    
  fNtuple->Branch("HitposAVRmmULadder",  HitposAVRmmULadder,  "HitposAVRmmULadder[HitNb]/F");
  fNtuple->Branch("HitposAVRmmVLadder",  HitposAVRmmVLadder,  "HitposAVRmmVLadder[HitNb]/F");
  fNtuple->Branch("HitmomMeVX",          HitmomMeVX,          "HitmomMeVX[HitNb]/F");
  fNtuple->Branch("HitmomMeVY",          HitmomMeVY,          "HitmomMeVY[HitNb]/F");
  fNtuple->Branch("HitmomMeVZ",          HitmomMeVZ,          "HitmomMeVZ[HitNb]/F");
  fNtuple->Branch("HitthetaLoc",         HitthetaLoc,         "HitthetaLoc[HitNb]/F");
  fNtuple->Branch("HitphiLoc",           HitphiLoc,           "HitphiLoc[HitNb]/F");
  fNtuple->Branch("HitglobalTime",       HitglobalTime,       "HitglobalTime[HitNb]/F");
  fNtuple->Branch("HitGeant4EdepMeV",    HitGeant4EdepMeV,    "HitGeant4EdepMeV[HitNb]/F");
  fNtuple->Branch("HitClusterSizeCol",   HitClusterSizeCol,   "HitClusterSizeCol[HitNb]/F");
  fNtuple->Branch("HitClusterSizeRow",   HitClusterSizeRow,   "HitClusterSizeRow[HitNb]/F");    
  fNtuple->Branch("HitRecoUmm",          HitRecoUmm,          "HitRecoUmm[HitNb]/F");
  fNtuple->Branch("HitRecoVmm",          HitRecoVmm,          "HitRecoVmm[HitNb]/F");    
  fNtuple->Branch("HitRecoULaddermm",    HitRecoULaddermm,    "HitRecoULaddermm[HitNb]/F");
  fNtuple->Branch("HitRecoVLaddermm",    HitRecoVLaddermm,    "HitRecoVLaddermm[HitNb]/F");
  fNtuple->Branch("HitPhiPrincipalAxis", HitPhiPrincipalAxis, "HitPhiPrincipalAxis[HitNb]/F");
  fNtuple->Branch("HitRMSPrincipalAxis", HitRMSPrincipalAxis, "HitRMSPrincipalAxis[HitNb]/F");
  fNtuple->Branch("HitThetaScattering",  HitThetaScattering,  "HitThetaScattering[HitNb]/F");
  fNtuple->Branch("HitNPixels",          HitNPixels,          "HitNPixels[HitNb]/I");
  fNtuple->Branch("HitFirstPixelIdx",    HitFirstPixelIdx,    "HitFirstPixelIdx[HitNb]/I");
    
  //Definition of Pixel Branch 
  fNtuple->Branch("PixelNb",          &PixelNb,           "PixelNb/I");
  fNtuple->Branch("PixelHitIdx",       PixelHitIdx,       "PixelHitIdx[PixelNb]/I");
  fNtuple->Branch("PixelGlobalIdx",    PixelGlobalIdx,    "PixelGlobalIdx[PixelNb]/I");
  fNtuple->Branch("PixelColumn",       PixelColumn,       "PixelColumn[PixelNb]/I");
  fNtuple->Branch("PixelRow",          PixelRow,          "PixelRow[PixelNb]w/I");
  fNtuple->Branch("PixelAnalogCharge", PixelAnalogCharge, "PixelAnalogCharge[PixelNb]/F");
  fNtuple->Branch("PixelSensorID",     PixelSensorID,     "PixelSensorID[PixelNb]/I");
  fNtuple->Branch("PixelUmm",          PixelUmm,          "PixelUmm[PixelNb]/F");
  fNtuple->Branch("PixelVmm",          PixelVmm,          "PixelVmm[PixelNb]/F");
  fNtuple->Branch("PixelStatus",       PixelStatus,       "PixelStatus[PixelNb]/I");

  //Definition of Saturation Branch    
  fNtuple->Branch("SaturationInfoNb",       &SatuNb,       "SaturationInfoNb/I");
  fNtuple->Branch("SaturationInfoLinIdx",    SatuLinIdx,   "SaturationInfoLinIdx[SaturationInfoNb]/I");
  fNtuple->Branch("SaturationInfoSensorID",  SatuSensorID, "SaturationInfoSensorID[SaturationInfoNb]/I");
   
  //Definition of EventSize Branch  
  fNtuple->Branch("EventSizeLadder1",  &EventSizeLadder1, "EventSizeLadder1/F");
  fNtuple->Branch("EventSizeLadder2",  &EventSizeLadder2, "EventSizeLadder2/F");

  if(FillNonSensitiveBranch) {
    //==============================
    //  Non-sensitive hits branches
    //==============================
    //Definition of Particles Branch
    fNtuple->Branch("NonSensitiveParticleNb",           &NonSensitiveParticleNb,             "NonSensitiveParticleNb/I");
    fNtuple->Branch("NonSensitiveParticleSensitiveIdx",  NonSensitiveParticleSensitiveIdx,   "NonSensitiveParticleSensitiveIdx[NonSensitiveParticleNb]/I");
    fNtuple->Branch("NonSensitiveParticlepdgID",         NonSensitiveParticlepdgID,          "NonSensitiveParticlepdgID[NonSensitiveParticleNb]/I");
    fNtuple->Branch("NonSensitiveParticleTrkVtX",        NonSensitiveParticleTrkVtX,         "NonSensitiveParticleTrkVtX[NonSensitiveParticleNb]/F");
    fNtuple->Branch("NonSensitiveParticleTrkVtY",        NonSensitiveParticleTrkVtY,         "NonSensitiveParticleTrkVtY[NonSensitiveParticleNb]/F");
    fNtuple->Branch("NonSensitiveParticleTrkVtZ",        NonSensitiveParticleTrkVtZ,         "NonSensitiveParticleTrkVtZ[NonSensitiveParticleNb]/F");
    fNtuple->Branch("NonSensitiveParticleNHits",         NonSensitiveParticleNHits,          "NonSensitiveParticleNHits[NonSensitiveParticleNb]/I");
    fNtuple->Branch("NonSensitiveParticleFirstHitIdx",   NonSensitiveParticleFirstHitIdx,    "NonSensitiveParticleFirstHitIdx[NonSensitiveParticleNb]/I");

    //Definition of Hits Branch
    fNtuple->Branch("NonSensitiveHitNb",                &NonSensitiveHitNb,                  "NonSensitiveHitNb/I");
    fNtuple->Branch("NonSensitiveHitParticleIdx",        NonSensitiveHitParticleIdx,         "NonSensitiveHitParticleIdx[NonSensitiveHitNb]/I");
    fNtuple->Branch("NonSensitiveHitposINmmX",           NonSensitiveHitposINmmX,            "NonSensitiveHitposINmmX[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitposINmmY",           NonSensitiveHitposINmmY,            "NonSensitiveHitposINmmY[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitposINmmZ",           NonSensitiveHitposINmmZ,            "NonSensitiveHitposINmmZ[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitglobalTimeINns",     NonSensitiveHitglobalTimeINns,      "NonSensitiveHitglobalTimeINns[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitmomINMeVX",          NonSensitiveHitmomINMeVX,           "NonSensitiveHitmomINMeVX[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitmomINMeVY",          NonSensitiveHitmomINMeVY,           "NonSensitiveHitmomINMeVY[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitmomINMeVZ",          NonSensitiveHitmomINMeVZ,           "NonSensitiveHitmomINMeVZ[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitEnergyINMeV",        NonSensitiveHitEnergyINMeV,         "NonSensitiveHitEnergyINMeV[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitposOUTmmX",          NonSensitiveHitposOUTmmX,           "NonSensitiveHitposOUTmmX[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitposOUTmmY",          NonSensitiveHitposOUTmmY,           "NonSensitiveHitposOUTmmY[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitposOUTmmZ",          NonSensitiveHitposOUTmmZ,           "NonSensitiveHitposOUTmmZ[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitmomOUTMeVX",         NonSensitiveHitmomOUTMeVX,          "NonSensitiveHitmomOUTMeVX[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitmomOUTMeVY",         NonSensitiveHitmomOUTMeVY,          "NonSensitiveHitmomOUTMeVY[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitmomOUTMeVZ",         NonSensitiveHitmomOUTMeVZ,          "NonSensitiveHitmomOUTMeVZ[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitEnergyOUTMeV",       NonSensitiveHitEnergyOUTMeV,        "NonSensitiveHitEnergyOUTMeV[NonSensitiveHitNb]/F");
    fNtuple->Branch("NonSensitiveHitGeant4EdepMeV",      NonSensitiveHitGeant4EdepMeV,       "NonSensitiveHitGeant4EdepMeV[NonSensitiveHitNb]/F");
  }

  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::save()
{ 
  
  // Saving the n-tuple to output file

  if(fRootFile) {
    //fRootFile->Write();      // Writing the histograms to the file
    fNtuple->Write();          // Writing the tree to the file
    fRootFile->Close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillMimosaArray(MimosaSimuTrackerHit* aHit, G4int index)
{

  // Fills up the hits in mimosa (sensitive volumes) list
  
  if(index == 0) HitsListInMimosa.clear();
  
  MyHit.pdgID         = aHit->GetTrackPDGCode();
  MyHit.sensorID      = aHit->GetMimosaSensorID() - 1;
  MyHit.trackID       = aHit->GetTrackID();
  MyHit.trackVtx      = (1.0/DistFactor)*aHit->GetTrackVertexPos();
  MyHit.energyDepMeV  = aHit->GetEdep()/MeV;
  
  MyHit.thetaLoc      = aHit->GetFirstStepLocalAngles().x()/AngleFactor;
  MyHit.phiLoc        = aHit->GetFirstStepLocalAngles().y()/AngleFactor;
  MyHit.posINmm_XYZ   = (1.0/DistFactor)*aHit->GetFirstStepPos();
  MyHit.posINmm_UVW   = (1.0/DistFactor)*aHit->GetFirstStepPosLocal();
  MyHit.momINGeV      = G4ThreeVector(aHit->GetFirstStep4Vector().px(),
				      aHit->GetFirstStep4Vector().py(),
				      aHit->GetFirstStep4Vector().pz());
  MyHit.momINGeV     *= (1.0/EPFactor);
  MyHit.EnergyINGeV   = (1.0/EPFactor)*aHit->GetFirstStep4Vector().e();
  MyHit.globalTimeIN  = aHit->GetFirstStepTime()/TimeFactor;
  MyHit.posOUTmm_XYZ  = (1.0/DistFactor)*aHit->GetLastStepPos();
  MyHit.posOUTmm_UVW  = (1.0/DistFactor)*aHit->GetLastStepPosLocal();
  MyHit.momOUTGeV     = G4ThreeVector(aHit->GetLastStep4Vector().px(),
				      aHit->GetLastStep4Vector().py(),
				      aHit->GetLastStep4Vector().pz());
  MyHit.momOUTGeV    *= (1.0/EPFactor);
  MyHit.EnergyOUTGeV  = (1.0/EPFactor)*aHit->GetLastStep4Vector().e();
  MyHit.globalTimeOUT = aHit->GetLastStepTime()/TimeFactor;
  
  HitsListInMimosa.push_back(MyHit);

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillNonSensitiveArray(MimosaSimuTrajectoryHit* aHit, G4int index)
{

  // Fills up the hits in non-sensitive volumes list
  
  if(!FillNonSensitiveBranch) return;
    
  if(index == 0) HitsListInNonSensitive.clear();
  
  MyHit.pdgID         = aHit->GetTrackPDGCode();
  MyHit.sensorID      = -1;
  MyHit.trackID       = aHit->GetTrackID();
  MyHit.trackVtx      = (1.0/DistFactor)*aHit->GetTrackVertexPos();
  MyHit.energyDepMeV  = aHit->GetEdep()/MeV;
  
  MyHit.thetaLoc      = -999.0;
  MyHit.phiLoc        = -999.0;
  MyHit.posINmm_XYZ   = (1.0/DistFactor)*aHit->GetFirstStepPos();
  MyHit.posINmm_UVW   = G4ThreeVector(-999.0,-999.0,-999.0);
  MyHit.momINGeV      = G4ThreeVector(aHit->GetFirstStep4Vector().px(),
				      aHit->GetFirstStep4Vector().py(),
				      aHit->GetFirstStep4Vector().pz());
  MyHit.momINGeV     *= (1.0/EPFactor);
  MyHit.EnergyINGeV   = (1.0/EPFactor)*aHit->GetFirstStep4Vector().e();
  MyHit.globalTimeIN  = aHit->GetFirstStepTime()/TimeFactor;
  MyHit.posOUTmm_XYZ  = (1.0/DistFactor)*aHit->GetLastStepPos();
  MyHit.posOUTmm_UVW  = G4ThreeVector(-999.0,-999.0,-999.0);
  MyHit.momOUTGeV     = G4ThreeVector(aHit->GetLastStep4Vector().px(),
				      aHit->GetLastStep4Vector().py(),
				      aHit->GetLastStep4Vector().pz());
  MyHit.momOUTGeV    *= (1.0/EPFactor);
  MyHit.EnergyOUTGeV  = (1.0/EPFactor)*aHit->GetLastStep4Vector().e();
  MyHit.globalTimeOUT = aHit->GetLastStepTime()/TimeFactor;
  
  HitsListInNonSensitive.push_back(MyHit);

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillTreeBlocks(void)
{

  EvtNumber++;
  
  // Fills up the n-tuple blocks
  
  //Fill list of particles producing hits in mimosas
  FillUpParticleList_Mimosa();
  
  //Fill list of particles producing hits in non-sensitive volumes
  FillUpParticleList_NonSensitive();
  
  HitsListInMimosa.clear();
  HitsListInNonSensitive.clear();

  // Digitize this event
  if(!NoDigitizationInAllSensitivePlanes) DigitizeEvent();

  // Fills hits in mimosa (sensitive volumes) block
  FillMimosaBlock();
  
  // Fills hits in non-sensitive volumes block
  if(FillNonSensitiveBranch) FillNonSensitiveBlock();

  // Fills line saturation block
  FillSaturationBlock();
  
  // Fills event size block
  FillEventSizeBlock();

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillUpParticleList_Mimosa(void)
{
 
  // Fill particle list out of hit list for Mimosa (sensitive volumes)
  // This function fills the particle list ListOfParticlesInMimosa from the list of all pixels HitsListInMimosa

  double AEdep_factor  = 1.0;
  AEdep_factor        *= Edep_correction_factor;
  
  //Get the list of particles in this event:
  ListOfParticlesInMimosa.clear();
  for(int ihit=0;ihit<int(HitsListInMimosa.size());ihit++) {
    int           trkID      = HitsListInMimosa[ihit].trackID;
    int           pdgID      = HitsListInMimosa[ihit].pdgID;
    G4ThreeVector trkVtx     = HitsListInMimosa[ihit].trackVtx;
    G4double      globalTime = HitsListInMimosa[ihit].globalTimeIN;

    //Check if particle is already included in the list
    bool IsNotAlreadyInList = true;
    for(int ipart=0;ipart<int(ListOfParticlesInMimosa.size());ipart++) {
      int trkID_tmp = ListOfParticlesInMimosa[ipart].trackID;
      
      if(trkID_tmp == trkID) {
	IsNotAlreadyInList = false;
	break;
      }
    }
    
    if(IsNotAlreadyInList) {
      AParticle_t aParticle;
      aParticle.pdgID      = pdgID;
      aParticle.trackID    = trkID;
      aParticle.trkVtx     = trkVtx;
      aParticle.globalTime = globalTime;
      aParticle.ListOfParticleHits.clear();
      ListOfParticlesInMimosa.push_back(aParticle);
    }
  }

  //Now collect the hits associated to a given particle in the list
  for(int ipart=0;ipart<int(ListOfParticlesInMimosa.size());ipart++) {
    int trkID_part = ListOfParticlesInMimosa[ipart].trackID;
    //int pdgID_part = ListOfParticlesInMimosa[ipart].pdgID;
    
    // Now loop over all the hits and fill only those produced by the current particle
    for(int ihit=0;ihit<int(HitsListInMimosa.size());ihit++) {
      int  trkID_hit = HitsListInMimosa[ihit].trackID;
     
      if(trkID_part != trkID_hit) continue;

      AParticleHits_t APartHit;
      APartHit.pdgID            = HitsListInMimosa[ihit].pdgID;
      APartHit.sensorID         = HitsListInMimosa[ihit].sensorID;
      APartHit.ladderID         = -1;
      APartHit.moduleID         = -1;
      APartHit.EdepMeV          = HitsListInMimosa[ihit].energyDepMeV*AEdep_factor;

      APartHit.posINmmXYZ       = HitsListInMimosa[ihit].posINmm_XYZ;
      APartHit.posINmmUVW       = HitsListInMimosa[ihit].posINmm_UVW;
      APartHit.momINGeV         = HitsListInMimosa[ihit].momINGeV;
      APartHit.EnergyINGeV      = HitsListInMimosa[ihit].EnergyINGeV;
      APartHit.thetaLoc         = HitsListInMimosa[ihit].thetaLoc;
      APartHit.phiLoc           = HitsListInMimosa[ihit].phiLoc;
      
      APartHit.posOUTmmXYZ      = HitsListInMimosa[ihit].posOUTmm_XYZ;
      APartHit.posOUTmmUVW      = HitsListInMimosa[ihit].posOUTmm_UVW;
      APartHit.momOUTGeV        = HitsListInMimosa[ihit].momOUTGeV;
      APartHit.EnergyOUTGeV     = HitsListInMimosa[ihit].EnergyOUTGeV;
      
      APartHit.globalTime       = HitsListInMimosa[ihit].globalTimeIN;
      
      APartHit.clusterSizeCol   = -1;
      APartHit.clusterSizeRow   = -1;
      APartHit.recoUVmm         = G4ThreeVector(-999.0,-999.0,-999.0);
      APartHit.recoUVLaddermm   = G4ThreeVector(-999.0,-999.0,-999.0);
      APartHit.phiPrincipalAxis = -999.0;
      APartHit.rmsPrincipalAxis = -999.0;
      APartHit.thetaScattering  = -999.0;
      
      APartHit.ListOfHitPixels.clear();
 
      ListOfParticlesInMimosa[ipart].ListOfParticleHits.push_back(APartHit);
    } // end of loop on hits

     //Order in time the hits of a particle:
     for(int iii=2;iii<=int(ListOfParticlesInMimosa[ipart].ListOfParticleHits.size());iii++) {
       for(int jjj=0;jjj<=int(ListOfParticlesInMimosa[ipart].ListOfParticleHits.size())-iii;jjj++) {
	 float time_jjj   = ListOfParticlesInMimosa[ipart].ListOfParticleHits[jjj].globalTime;
	 float time_jjjp1 = ListOfParticlesInMimosa[ipart].ListOfParticleHits[jjj+1].globalTime;
	 
	 if(time_jjj > time_jjjp1) {
	   AParticleHits_t APartHit_aux;
	   EqualizeAParticleHits(APartHit_aux,                                             ListOfParticlesInMimosa[ipart].ListOfParticleHits[jjj]);
	   EqualizeAParticleHits(ListOfParticlesInMimosa[ipart].ListOfParticleHits[jjj],   ListOfParticlesInMimosa[ipart].ListOfParticleHits[jjj+1]);
	   EqualizeAParticleHits(ListOfParticlesInMimosa[ipart].ListOfParticleHits[jjj+1], APartHit_aux);
	}
      }
    }
    
    ListOfParticlesInMimosa[ipart].trkVtx     = ListOfParticlesInMimosa[ipart].ListOfParticleHits[0].posINmmXYZ;
    ListOfParticlesInMimosa[ipart].globalTime = ListOfParticlesInMimosa[ipart].ListOfParticleHits[0].globalTime;
    
    if(!FillNonSensitiveBranch) continue;
    
    for(int ihit=0;ihit<int(HitsListInNonSensitive.size());ihit++) {
      int  trkID_hit = HitsListInNonSensitive[ihit].trackID;
      
      if(trkID_part != trkID_hit) continue;
      
      if(ListOfParticlesInMimosa[ipart].globalTime > HitsListInNonSensitive[ihit].globalTimeIN) {
	ListOfParticlesInMimosa[ipart].trkVtx     = HitsListInNonSensitive[ihit].posINmm_XYZ;
	ListOfParticlesInMimosa[ipart].globalTime = HitsListInNonSensitive[ihit].globalTimeIN;
      }
      
    }
    
  } // End of particle hits loop

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillUpParticleList_NonSensitive(void)
{  

  // Fill particle list out of hit list for non-sensitive volumes
  // This function fills the particle list ListOfParticlesInNonSensitive from the list of all pixels HitsListInNonSensitive

  if(!FillNonSensitiveBranch) return;
  
  //Get the list of particles in this event:
  ListOfParticlesInNonSensitive.clear();
  for(int ihit=0;ihit<int(HitsListInNonSensitive.size());ihit++) {
    int           trkID      = HitsListInNonSensitive[ihit].trackID;
    int           pdgID      = HitsListInNonSensitive[ihit].pdgID;
    G4ThreeVector trkVtx     = HitsListInNonSensitive[ihit].trackVtx;
    G4double      globalTime = HitsListInNonSensitive[ihit].globalTimeIN;

    //Check if particle is already included in the list
    bool IsNotAlreadyInList = true;
    for(int ipart=0;ipart<int(ListOfParticlesInNonSensitive.size());ipart++) {
      int trkID_tmp = ListOfParticlesInNonSensitive[ipart].trackID;
      
      if(trkID_tmp == trkID) {
	IsNotAlreadyInList = false;
	break;
      }
    }
    
    if(IsNotAlreadyInList) {
      AParticle_t aParticle;
      aParticle.pdgID      = pdgID;
      aParticle.trackID    = trkID;
      aParticle.trkVtx     = trkVtx;
      aParticle.globalTime = globalTime;
      aParticle.ListOfParticleHits.clear();
      ListOfParticlesInNonSensitive.push_back(aParticle);
    }
  }

  //Now collect the hits associated to a given particle in the list
  for(int ipart=0;ipart<int(ListOfParticlesInNonSensitive.size());ipart++) {
    int trkID_part = ListOfParticlesInNonSensitive[ipart].trackID;
    //int pdgID_part = ListOfParticlesInNonSensitive[ipart].pdgID;
    
    // Now loop over all the hits and fill only those produced by the current particle
    for(int ihit=0;ihit<int(HitsListInNonSensitive.size());ihit++) {
      int  trkID_hit = HitsListInNonSensitive[ihit].trackID;
     
      if(trkID_part != trkID_hit) continue;

      AParticleHits_t APartHit;
      APartHit.pdgID            = HitsListInNonSensitive[ihit].pdgID;
      APartHit.sensorID         = HitsListInNonSensitive[ihit].sensorID;
      APartHit.ladderID         = -1;
      APartHit.moduleID         = -1;
      APartHit.EdepMeV          = HitsListInNonSensitive[ihit].energyDepMeV;

      APartHit.posINmmXYZ       = HitsListInNonSensitive[ihit].posINmm_XYZ;
      APartHit.posINmmUVW       = HitsListInNonSensitive[ihit].posINmm_UVW;
      APartHit.momINGeV         = HitsListInNonSensitive[ihit].momINGeV;
      APartHit.EnergyINGeV      = HitsListInNonSensitive[ihit].EnergyINGeV;
      APartHit.thetaLoc         = HitsListInNonSensitive[ihit].thetaLoc;
      APartHit.phiLoc           = HitsListInNonSensitive[ihit].phiLoc;
      
      APartHit.posOUTmmXYZ      = HitsListInNonSensitive[ihit].posOUTmm_XYZ;
      APartHit.posOUTmmUVW      = HitsListInNonSensitive[ihit].posOUTmm_UVW;
      APartHit.momOUTGeV        = HitsListInNonSensitive[ihit].momOUTGeV;
      APartHit.EnergyOUTGeV     = HitsListInNonSensitive[ihit].EnergyOUTGeV;
      
      APartHit.globalTime       = HitsListInNonSensitive[ihit].globalTimeIN;
      
      APartHit.clusterSizeCol   = -1;
      APartHit.clusterSizeRow   = -1;
      APartHit.recoUVmm         = G4ThreeVector(-999.0,-999.0,-999.0);
      APartHit.recoUVLaddermm   = G4ThreeVector(-999.0,-999.0,-999.0);
      APartHit.phiPrincipalAxis = -999.0;
      APartHit.rmsPrincipalAxis = -999.0;
      APartHit.thetaScattering  = -999.0;
      
      APartHit.ListOfHitPixels.clear();
 
      ListOfParticlesInNonSensitive[ipart].ListOfParticleHits.push_back(APartHit);
    } // end of loop on hits

     //Order in time the hits of a particle:
     for(int iii=2;iii<=int(ListOfParticlesInNonSensitive[ipart].ListOfParticleHits.size());iii++) {
       for(int jjj=0;jjj<=int(ListOfParticlesInNonSensitive[ipart].ListOfParticleHits.size())-iii;jjj++) {
	 float time_jjj   = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[jjj].globalTime;
	 float time_jjjp1 = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[jjj+1].globalTime;
	 
	 if(time_jjj > time_jjjp1) {
	   AParticleHits_t APartHit_aux;
	   EqualizeAParticleHits(APartHit_aux,                                             ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[jjj]);
	   EqualizeAParticleHits(ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[jjj],   ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[jjj+1]);
	   EqualizeAParticleHits(ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[jjj+1], APartHit_aux);
	}
      }
    }
    
    ListOfParticlesInNonSensitive[ipart].trkVtx     = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[0].posINmmXYZ;
    ListOfParticlesInNonSensitive[ipart].globalTime = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[0].globalTime;
    
    for(int ihit=0;ihit<int(HitsListInMimosa.size());ihit++) {
      int  trkID_hit = HitsListInMimosa[ihit].trackID;
      
      if(trkID_part != trkID_hit) continue;
      
      if(ListOfParticlesInNonSensitive[ipart].globalTime > HitsListInMimosa[ihit].globalTimeIN) {
	ListOfParticlesInNonSensitive[ipart].trkVtx     = HitsListInMimosa[ihit].posINmm_XYZ;
	ListOfParticlesInNonSensitive[ipart].globalTime = HitsListInMimosa[ihit].globalTimeIN;
      }
      
    }
    
  } // End of particle hits loop
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillUpParticleList(std::vector<Hit_t>  HitsList, std::vector<AParticle_t>&  ListOfParticles)
{

  // This function fills the particle list ListOfParticles from the list of all pixels HitsList

  //Get the list of particles in this event:
  ListOfParticles.clear();
  for(int ihit=0;ihit<int(HitsList.size());ihit++) {
    int           trkID   = HitsList[ihit].trackID;
    int           pdgID   = HitsList[ihit].pdgID;
    G4ThreeVector trkVtx  = HitsList[ihit].trackVtx;

    //Check if particle is already included in the list
    bool IsNotAlreadyInList = true;
    for(int ipart=0;ipart<int(ListOfParticles.size());ipart++) {
      int trkID_tmp = ListOfParticles[ipart].trackID;
      
      if(trkID_tmp == trkID) {
	IsNotAlreadyInList = false;
	break;
      }
    }
    
    if(IsNotAlreadyInList) {
      AParticle_t aParticle;
      aParticle.pdgID   = pdgID;
      aParticle.trackID = trkID;
      aParticle.trkVtx  = trkVtx;
      aParticle.ListOfParticleHits.clear();
      ListOfParticles.push_back(aParticle);
    }
  }

  //Now collect the hits associated to a given particle in the list
  for(int ipart=0;ipart<int(ListOfParticles.size());ipart++) {
    int trkID_part = ListOfParticles[ipart].trackID;
    //int pdgID_part = ListOfParticles[ipart].pdgID;
    
    // Now loop over all the hits and fill only those produced by the current particle
    for(int ihit=0;ihit<int(HitsList.size());ihit++) {
      int  trkID_hit = HitsList[ihit].trackID;
      
      if(trkID_part != trkID_hit) continue;

      AParticleHits_t APartHit;
      APartHit.pdgID            = HitsList[ihit].pdgID;
      APartHit.sensorID         = HitsList[ihit].sensorID;
      APartHit.ladderID         = -1;
      APartHit.moduleID         = -1;
      APartHit.EdepMeV          = HitsList[ihit].energyDepMeV;

      APartHit.posINmmXYZ       = HitsList[ihit].posINmm_XYZ;
      APartHit.posINmmUVW       = HitsList[ihit].posINmm_UVW;
      APartHit.momINGeV         = HitsList[ihit].momINGeV;
      APartHit.EnergyINGeV      = HitsList[ihit].EnergyINGeV;
      APartHit.thetaLoc         = HitsList[ihit].thetaLoc;
      APartHit.phiLoc           = HitsList[ihit].phiLoc;
      
      APartHit.posOUTmmXYZ      = HitsList[ihit].posOUTmm_XYZ;
      APartHit.posOUTmmUVW      = HitsList[ihit].posOUTmm_UVW;
      APartHit.momOUTGeV        = HitsList[ihit].momOUTGeV;
      APartHit.EnergyOUTGeV     = HitsList[ihit].EnergyOUTGeV;
      
      APartHit.globalTime       = HitsList[ihit].globalTimeIN;
      
      APartHit.clusterSizeCol   = -1;
      APartHit.clusterSizeRow   = -1;
      APartHit.recoUVmm         = G4ThreeVector(-999.0,-999.0,-999.0);
      APartHit.recoUVLaddermm   = G4ThreeVector(-999.0,-999.0,-999.0);
      APartHit.phiPrincipalAxis = -999.0;
      APartHit.rmsPrincipalAxis = -999.0;
      APartHit.thetaScattering  = -999.0;
      
      APartHit.ListOfHitPixels.clear();
 
      ListOfParticles[ipart].ListOfParticleHits.push_back(APartHit);

      //if(trkID_part == trkID_hit) break;
    } // end of loop on hits

     //Order in time the hits of a particle:
     for(int iii=2;iii<=int(ListOfParticles[ipart].ListOfParticleHits.size());iii++) {
       for(int jjj=0;jjj<=int(ListOfParticles[ipart].ListOfParticleHits.size())-iii;jjj++) {
	 float time_jjj   = ListOfParticles[ipart].ListOfParticleHits[jjj].globalTime;
	 float time_jjjp1 = ListOfParticles[ipart].ListOfParticleHits[jjj+1].globalTime;
	 
	 if(time_jjj > time_jjjp1) {
	   AParticleHits_t APartHit_aux;
	   EqualizeAParticleHits(APartHit_aux,                                     ListOfParticles[ipart].ListOfParticleHits[jjj]);
	   EqualizeAParticleHits(ListOfParticles[ipart].ListOfParticleHits[jjj],   ListOfParticles[ipart].ListOfParticleHits[jjj+1]);
	   EqualizeAParticleHits(ListOfParticles[ipart].ListOfParticleHits[jjj+1], APartHit_aux);
	}
      }
    }

    
  } // End of particle hits loop
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::DigitizeEvent(void)
{

  // This function digitizes the event. loop over all the hits in mimosa (sensitive volumes)
  // and digitizes each one. For each plane will use the digitization configuration specified 
  // in the corresponding plane block in the config file
  
  // Loops over all the hist in mimosa
  for(int ipart=0;ipart<int(ListOfParticlesInMimosa.size());ipart++) {
    for(int ihit=0;ihit<int(ListOfParticlesInMimosa[ipart].ListOfParticleHits.size());ihit++) {
      
      // For a given hit check the index of this sensitive plane in the ListOfSensitivePlanes
      int idx_SD = GetIdxInSensitivePlaneList(ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].sensorID);
      if(idx_SD == -1) continue;

      // If digitization if turned on for this plane digitize the hit
      if(ListOfSensitivePlanes[idx_SD].digit->GetDoDigitization()) DigitizeHit(idx_SD, ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit]);
      
    }
  }
  
  // Once all the hits have been digitized, proceed t producing pixels from noise
  FillPixelsFromNoise(ListOfPixelsFromNoise);

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::DigitizeHit(int idx_SD, AParticleHits_t& aHit)
{
  
  // This function digitizes a hit using the specification in ListOfSensitivePlanes[idx_SD]
 
  //Check if particle is charged. If not, do nothing.
  if(!IsParticleCharged(aHit.pdgID)) return;
  
  // Threshold in electrons units
  Double_t Discri_Threshold_electrons = ListOfSensitivePlanes[idx_SD].digit->GetNoiseElectrons()*ListOfSensitivePlanes[idx_SD].DigitizeThre;
  
  // Segment size for charge transport as a fraction of epitaxial layer thickness
  Double_t mm_to_um = 1.0e+3; // comvertion factor from mm to um
  Double_t InputUpos  = aHit.posINmmUVW.x()  * mm_to_um;
  Double_t InputVpos  = aHit.posINmmUVW.y()  * mm_to_um;
  Double_t OutputUpos = aHit.posOUTmmUVW.x() * mm_to_um;
  Double_t OutputVpos = aHit.posOUTmmUVW.y() * mm_to_um;
  
  Double_t epsilon_epi = 1.0e-2;
  Double_t TheCurrentEpi = TMath::Abs(aHit.posINmmUVW.z() - aHit.posOUTmmUVW.z()) * mm_to_um;
  Double_t TheNominalEpi = ListOfSensitivePlanes[idx_SD].digit->GetEffectiveEpi();
  if(TheCurrentEpi < epsilon_epi) {
    TheCurrentEpi = TheNominalEpi;
    cout << "WARNING: The current epi thickness is zero, setting it to the nominal epi thickness (" << TheNominalEpi << " um). " << endl;
    cout << "         (W-in,W-out) = (" << aHit.posINmmUVW.z() * mm_to_um << "," << aHit.posOUTmmUVW.z() * mm_to_um << ") um " << endl;
  }
  if(TMath::Abs(TheCurrentEpi - TheNominalEpi) > epsilon_epi) {
    cout << "WARNING: The current epi thickness (" << TheCurrentEpi << " um) is different from the nominal epi thickness (" << TheNominalEpi << " um). " << endl;
    cout << "         (W-in,W-out) = (" << aHit.posINmmUVW.z() * mm_to_um << "," << aHit.posOUTmmUVW.z() * mm_to_um << ") um " << endl;
    ListOfSensitivePlanes[idx_SD].digit->SetEffectiveEpi(TheCurrentEpi);
    cout << "         The new epi-thickness is = " << ListOfSensitivePlanes[idx_SD].digit->GetEffectiveEpi() << " um " << endl;
  }
  Double_t SegmentSize  = ListOfSensitivePlanes[idx_SD].digit->GetEffectiveEpi()/20.0;

  double TheDepositedEGeant4  = aHit.EdepMeV*1.0e+6/3.6;   // deposited energy by ionization in electrons
  //TheDepositedEGeant4        *= Edep_correction_factor;

  //If TheDepositedEGeant4 negative the digitization code will generate its own deposited energy
  //using a Landau distribution with MPV = 80 e-/um and Width = 18 e-/um
  //TheDepositedEGeant4 = -1;

  // perform the deposited charge transportation and digitization
  std::vector<int>   fPixelMap;          // index of pixels which collects some charge
  std::vector<float> fAnalogChargeMap;   // Analog charge  in e- units
  std::vector<int>   fDigitalChargeMap;  // digital charge in ADC units
  ListOfSensitivePlanes[idx_SD].digit->Transport_Noise_Digitisation(SegmentSize,
                                                                    InputUpos,  InputVpos,
								    OutputUpos, OutputVpos,
								    TheDepositedEGeant4,
								    aHit.thetaLoc,aHit.phiLoc,
								    fPixelMap,fAnalogChargeMap,fDigitalChargeMap,
								    Discri_Threshold_electrons);
  if(TMath::Abs(TheCurrentEpi - TheNominalEpi) > epsilon_epi) {
    ListOfSensitivePlanes[idx_SD].digit->SetEffectiveEpi(TheNominalEpi);
    cout << "         Restoring the nominal epi-thickness is = " << ListOfSensitivePlanes[idx_SD].digit->GetEffectiveEpi() << " um " << endl;
  }

  // Now fill up this hit pixel list
  aHit.ListOfHitPixels.clear();
  for(int ipix=0;ipix<int(fPixelMap.size());ipix++) {
    // if digital output sensor, then will only fill in list pixels which passed threshold
    if(ListOfSensitivePlanes[idx_SD].digit->GetAnalysisMode() == 3 && fDigitalChargeMap[ipix] == 0) continue;
    
    int Row,Col;
    if(fPixelMap[ipix] < 0) cout << "WARNING: pixel map is negative = " << fPixelMap[ipix] << endl;
    ListOfSensitivePlanes[idx_SD].digit->GetXYPixelNumber(Col, Row, fPixelMap[ipix]);
    
    if(Col < 0 || Col > ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU() - 1) {
      cout << "  WARNING: column " << Col << " of pixel from hit is outside range (0," << ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU()-1 << ")" << endl;
      continue;
    }
    if(Row < 0 || Row > ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV() - 1) {
      cout << "  WARNING: column " << Row << " of pixel from hit is outside range (0," << ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV()-1 << ")" << endl;
      continue;
    }
   
    double U = 0.0;
    double V = 0.0;
    ListOfSensitivePlanes[idx_SD].digit->ComputePixelPositionUV_FromColRow(Col, Row, U, V);
    U /= mm_to_um;
    V /= mm_to_um;

    double TheAnalogCherge  = fAnalogChargeMap[ipix];                                               //analog charge in electrons
    TheAnalogCherge        *= ListOfSensitivePlanes[idx_SD].digit->GetCalib();                      // Use calibration to conver noise from electrons to volts
    TheAnalogCherge         = ListOfSensitivePlanes[idx_SD].digit->ConvertToADC(TheAnalogCherge);   // Convert noise in volts to ADC units
    
    AHitPixels_t ThePixel;
    ThePixel.globalIdx       = fPixelMap[ipix];
    ThePixel.col             = Col;
    ThePixel.row             = Row;
    ThePixel.ChargeAnalog    = TheAnalogCherge;
    ThePixel.sensorID        = aHit.sensorID;
    ThePixel.PosUVmm         = G4ThreeVector(U,V,0.0);
    ThePixel.status          = 1;

    aHit.ListOfHitPixels.push_back(ThePixel);
  } //end of pixel loop

  return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillPixelsFromNoise(std::vector<AHitPixels_t> &ListOfPixels)
{

  ListOfPixelsFromNoise.clear();
  
  for(int iplane=0;iplane<int(ListOfSensitivePlanes.size());iplane++) {
    if(ListOfSensitivePlanes[iplane].digit->GetAnalysisMode() == 2)       FillPixelsFromNoise_Analog(iplane,ListOfPixels);
    else if(ListOfSensitivePlanes[iplane].digit->GetAnalysisMode() == 3)  FillPixelsFromNoise_Digital(iplane,ListOfPixels);
  }

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillPixelsFromNoise_Analog(int idx_SD, std::vector<AHitPixels_t> &ListOfPixels)
{

  // This function fills up the list of pixels from noise for Analog output sensors. 
  // It Will use the Noise of the sensor (in ADC units) to fill up the pixels of this 
  // sensor which where not turned on by a particle
  
  // Conver noise in electons to ADC units
  Double_t TheNoise  = ListOfSensitivePlanes[idx_SD].digit->GetNoiseElectrons();      // noise in electons
  TheNoise          /= ListOfSensitivePlanes[idx_SD].digit->GetCalib();               // Use calibration to conver noise from electrons to volts
  TheNoise           = ListOfSensitivePlanes[idx_SD].digit->ConvertToADC(TheNoise);   // Convert noise in volts to ADC units
  
  //Loop over all the pixels of this sensor
  for(int Col=0;Col<ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU();Col++) {
    for(int Row=0;Row<ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV();Row++) {

      if(Col < 0 || Col > ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU() - 1) {
        cout << "  WARNING: column " << Col << " of pixel from hit is outside range (0," << ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU()-1 << ")" << endl;
        continue;
      }
      if(Row < 0 || Row > ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV() - 1) {
        cout << "  WARNING: column " << Row << " of pixel from hit is outside range (0," << ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV()-1 << ")" << endl;
        continue;
      }
      
      //Check if a pixel with this coordinates (Col,Row) and for this sensor is already in the pixel list
      if(IsAlreadyInPixelList(Col,Row,ListOfSensitivePlanes[idx_SD].sensorID)) continue;
      
      // Fill the pixel from noise list
      
      double U = 0;
      double V = 0;
      ListOfSensitivePlanes[idx_SD].digit->ComputePixelPositionUV_FromColRow(Col, Row, U, V);
      U /= 1.0e+3;
      V /= 1.0e+3;
      
      AHitPixels_t ThePixelNoise;
      ThePixelNoise.globalIdx             = Col + ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU()*Row;
      ThePixelNoise.col                   = Col;
      ThePixelNoise.row                   = Row;
      ThePixelNoise.ChargeAnalog          = rand->Gaus(0.0,TheNoise);
      ThePixelNoise.sensorID              = ListOfSensitivePlanes[idx_SD].sensorID;
      ThePixelNoise.PosUVmm               = G4ThreeVector(U,V,0);
      ThePixelNoise.status                = 1;
      
      ListOfPixels.push_back(ThePixelNoise);
    }  // end of loop over columns 
  }  // end of loop over rows

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillPixelsFromNoise_Digital(int idx_SD, std::vector<AHitPixels_t> &ListOfPixels)
{

  // This function fills up the list of pixels from noise for digitial output sensors. 
  // It Will use the Noise of occupancy variable in the config file to estimate the number of pixels 
  // turned on by noise in average and generate by Poisson distribution the actual number of fake hits.
  // Will then position them with an uniform distribution in the matrix surface.

  // Get average number of fake hits
  double NPixAvrNoise = ListOfSensitivePlanes[idx_SD].DigitizeOcc*ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU()*ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV();
  // Generate the number of fake hits with a Poisson distribution
  int NPixelNoise     = rand->Poisson(NPixAvrNoise);

  
  // Now loop over the fake hits and place then randomly with a uniform distribution in pixel matrix surface. 
  for(int ipix=0;ipix<NPixelNoise;ipix++){
    int Col  = RoundOffInt(rand->Uniform(-0.5,ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU() - 0.5));
    int Row  = RoundOffInt(rand->Uniform(-0.5,ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU() - 0.5));
    
    if(Col < 0 || Col > ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU() - 1) {
      cout << "  WARNING: column " << Col << " of pixel from hit is outside range (0," << ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU()-1 << ")" << endl;
      continue;
    }
    if(Row < 0 || Row > ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV() - 1) {
      cout << "  WARNING: column " << Row << " of pixel from hit is outside range (0," << ListOfSensitivePlanes[idx_SD].digit->GetNpixelsV()-1 << ")" << endl;
      continue;
    }
    
    //Check if a pixel with this coordinates (Col,Row) and for this sensor is already in the pixel list
    if(IsAlreadyInPixelList(Col,Row,ListOfSensitivePlanes[idx_SD].sensorID)) continue;

    double U = 0;
    double V = 0;
    ListOfSensitivePlanes[idx_SD].digit->ComputePixelPositionUV_FromColRow(Col, Row, U, V);
    U /= 1.0e+3;
    V /= 1.0e+3;

    AHitPixels_t ThePixelNoise;
    ThePixelNoise.globalIdx             = Col + ListOfSensitivePlanes[idx_SD].digit->GetNpixelsU()*Row;
    ThePixelNoise.col                   = Col;
    ThePixelNoise.row                   = Row;
    ThePixelNoise.ChargeAnalog          = -10000;
    ThePixelNoise.sensorID              = ListOfSensitivePlanes[idx_SD].sensorID;
    ThePixelNoise.PosUVmm               = G4ThreeVector(U,V,0);
    ThePixelNoise.status                = 1;
       
    ListOfPixels.push_back(ThePixelNoise);
  }

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::InitSensitivePlaneList(void)
{

  // This function initializes and fill the list of sensitive planes. 
  // Also initializes the digitizer object
  
  NoDigitizationInAllSensitivePlanes = true;
  //if(NoDigitizationInAllSensitivePlanes) cout << "NoDigitizationInAllSensitivePlanes is true" << endl; 
  
  ListOfSensitivePlanes.clear();
  
  cout << endl;
  cout << "Number of sensitive planes = " << fSetup->GetNSensitivePlanes() << endl;
  for(int iplane=0;iplane<fSetup->GetNSensitivePlanes();iplane++) {
    ASensitivePlane_t APlaneSD;
    APlaneSD.sensorID       = fSetup->GetSensitivePlanesIdx(iplane) - 1;
    APlaneSD.DigitizeOcc    = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeOcc;
    APlaneSD.DigitizeThre   = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeThre;

    double aPitchU          = fSetup->GetPlanePar(APlaneSD.sensorID+1).Pitch[0];
    double aPitchV          = fSetup->GetPlanePar(APlaneSD.sensorID+1).Pitch[1];
    int    aNpixelsU        = fSetup->GetPlanePar(APlaneSD.sensorID+1).Strips[0];
    int    aNpixelsV        = fSetup->GetPlanePar(APlaneSD.sensorID+1).Strips[1];
    double aEffectiveEpi    = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneThickness;
    aEffectiveEpi          *= fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneEpiThickness;
    double ashiftU          = 0.0;
    double ashiftV          = 0.0;
    double aNoiseElectrons  = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeNoise;
    int    aMapping         = fSetup->GetPlanePar(APlaneSD.sensorID+1).Mapping;
    int    aAnalysisMode    = fSetup->GetPlanePar(APlaneSD.sensorID+1).AnalysisMode;
    TString aTransportModel = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitization;
    double aCalib           = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeCalib;
    int    aADCbits         = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeADCbits;
    double aADCRangeMin     = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeADCMin;
    double aADCRangeMax     = fSetup->GetPlanePar(APlaneSD.sensorID+1).PlaneDigitizeADCMax;

    cout << " - Sensitive plane " << iplane+1 << " has ID = " << fSetup->GetSensitivePlanesIdx(iplane) << endl;
    cout << "   * Mapping        = " << aMapping << endl;
    cout << "   * AnalysisMode   = " << aAnalysisMode << endl;
    cout << "   * PitchU/V       = " << aPitchU   << "/" << aPitchV   << " um" << endl;
    cout << "   * NpixelsU/V     = " << aNpixelsU << "/" << aNpixelsV << endl;
    cout << "   * Epilayer       = " << aEffectiveEpi << " um" << endl;
    cout << "   * Noise          = " << aNoiseElectrons << " elec" << endl;
    cout << "   * TransportMocel = " << aTransportModel.Data() << endl;
    if(aAnalysisMode == 3) {
      cout << "   * Occupancy      = " << APlaneSD.DigitizeOcc  << endl;
      cout << "   * Threshold      = " << APlaneSD.DigitizeThre << endl;
    }
    else if(aAnalysisMode == 2) {
      cout << "   * Calib          = " << aCalib << " elec/Volt" << endl;
      cout << "   * ADCbits        = " << aADCbits << endl;
      cout << "   * ADC Range      = (" << aADCRangeMin << "," << aADCRangeMax << ") Volts" << endl;
    }
    
    APlaneSD.digit = new DIGIT_B2();
    APlaneSD.digit->SetSensorParameters(aPitchU,   aPitchV,
				        ashiftU,   ashiftV,
				        aNpixelsU, aNpixelsV,
				        aEffectiveEpi, 
				        aNoiseElectrons,
				        aMapping,
				        aAnalysisMode,
				        aTransportModel,
				        aCalib,
				        aADCbits,aADCRangeMin,aADCRangeMax);
    APlaneSD.digit->SetGlobalSeed(fSetup->GetAnalysisPar().MCSeed);
    if(APlaneSD.digit->GetDoDigitization()) cout << "   * Will perform digitization for this plane." << endl;
    else                                    cout << "   * Will not perform digitization for this plane." << endl;
    cout << endl;
    if(NoDigitizationInAllSensitivePlanes && APlaneSD.digit->GetDoDigitization()) {
      //cout << "NoDigitizationInAllSensitivePlanes is true && Plane digitization is true. Changed to false" << endl;
      NoDigitizationInAllSensitivePlanes = false;
    }
    //if(NoDigitizationInAllSensitivePlanes) cout << "NoDigitizationInAllSensitivePlanes is true" << endl; 
    
    ListOfSensitivePlanes.push_back(APlaneSD);
  }
  cout << endl;
  
  if(NoDigitizationInAllSensitivePlanes) {
    cout << "NO DIGITIZATION WILL BE PERFORMED FOR ANY PLANE!!!" << endl;
    cout << endl;
  }
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillMimosaBlock(void)
{

  // This function fill the Mimosa (sensitive volumes) hits block 

  // Filling up the particles with their hits and their turned on pixels
  ParticleNb       = int(ListOfParticlesInMimosa.size());
  if(ParticleNb > NMaxParticle) {
    cout << endl;
    cout << "ParticleNb variable with value " << ParticleNb << " goes beyong the maximum allowed value " << NMaxParticle << ". Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  int HitNbCount   = 0;
  int PixelNbCount = 0;
  for(int ipart=0;ipart<int(ListOfParticlesInMimosa.size());ipart++) { 
    ParticleBKGType[ipart] = -1;
    ParticletrkID[ipart]   = ListOfParticlesInMimosa[ipart].trackID;
    ParticlepdgID[ipart]   = ListOfParticlesInMimosa[ipart].pdgID;
    ParticleTrkVtX[ipart]  = ListOfParticlesInMimosa[ipart].trkVtx.x();
    ParticleTrkVtY[ipart]  = ListOfParticlesInMimosa[ipart].trkVtx.y();
    ParticleTrkVtZ[ipart]  = ListOfParticlesInMimosa[ipart].trkVtx.z();
    ParticleNHits[ipart]   = ListOfParticlesInMimosa[ipart].ListOfParticleHits.size();
    if(ParticleNHits[ipart] == 0) ParticleFirstHitIdx[ipart] = -1;
    else                          ParticleFirstHitIdx[ipart] = HitNbCount;

    for(int ihit=0;ihit<int(ListOfParticlesInMimosa[ipart].ListOfParticleHits.size());ihit++) {
      if(HitNbCount+1 > NMaxHit) {
        cout << endl;
        cout << "Number of hits in Mimosa with value " << HitNbCount+1 << " goes beyong the maximum allowed value " << NMaxHit << ". Exiting now!!!" << endl;
        cout << endl;
        assert(false);
      }

      HitParticleIdx[HitNbCount]      = ipart;
      HitsensorID[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].sensorID;
      HitladderID[HitNbCount]         = -1;
      HitmoduleID[HitNbCount]         = -1;

      HitposINmmX[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmXYZ.x();
      HitposINmmY[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmXYZ.y();
      HitposINmmZ[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmXYZ.z();
      
      HitposINmmU[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmUVW.x();
      HitposINmmV[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmUVW.y();
      HitposINmmW[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmUVW.z();
      
      HitposOUTmmU[HitNbCount]        = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posOUTmmUVW.x();
      HitposOUTmmV[HitNbCount]        = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posOUTmmUVW.y();
      HitposOUTmmW[HitNbCount]        = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posOUTmmUVW.z();
      
      HitposAVRmmU[HitNbCount]        = (ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmUVW.x() + ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posOUTmmUVW.x())*0.5;
      HitposAVRmmV[HitNbCount]        = (ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmUVW.y() + ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posOUTmmUVW.y())*0.5;
      HitposAVRmmW[HitNbCount]        = TMath::Abs(ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posINmmUVW.z() - ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].posOUTmmUVW.z());

      HitposAVRmmULadder[HitNbCount]  = -999.0;
      HitposAVRmmVLadder[HitNbCount]  = -999.0;
      
      HitmomMeVX[HitNbCount]          = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].momINGeV.x()*1.0e+3;
      HitmomMeVY[HitNbCount]          = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].momINGeV.y()*1.0e+3;
      HitmomMeVZ[HitNbCount]          = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].momINGeV.z()*1.0e+3;

      HitthetaLoc[HitNbCount]         = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].thetaLoc;
      HitphiLoc[HitNbCount]           = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].phiLoc;
      HitglobalTime[HitNbCount]       = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].globalTime;
      HitGeant4EdepMeV[HitNbCount]    = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].EdepMeV;

      HitNPixels[HitNbCount]          = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels.size();
      if(HitNPixels[HitNbCount]==0) HitFirstPixelIdx[HitNbCount] = -1;
      else 			    HitFirstPixelIdx[HitNbCount] = PixelNbCount;
      
      HitClusterSizeCol[HitNbCount]   = -1;
      HitClusterSizeRow[HitNbCount]   = -1;
      HitRecoUmm[HitNbCount]          = -1.0e+5;
      HitRecoVmm[HitNbCount]          = -1.0e+5;
      HitRMSPrincipalAxis[HitNbCount] = -1; 
      HitPhiPrincipalAxis[HitNbCount] = -1000;
      HitThetaScattering[HitNbCount]  = -999.0;
      
      int idx_SD = GetIdxInSensitivePlaneList(ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].sensorID);
      if(idx_SD == -1) cout << "  WARNING, Inside FillMimosaBlock functon. The index of the sensitive plane is -1!!!!!" << endl;
      
      // If digitization is not turned on the rest of the variables are not filled and dummy values are asigned
      if(!ListOfSensitivePlanes[idx_SD].digit->GetDoDigitization()) HitNPixels[HitNbCount] = -1;
      else {
        int Colmax = -100000;
        int Colmin = +100000;
        int Rowmax = -100000;
        int Rowmin = +100000;
      
        float SumPosU  = 0;
        float SumPosV  = 0;
          
        float SumPosUU = 0;
        float SumPosVV = 0;
        float SumPosUV = 0;

        float SumOfWeights = 0.0;
        int MyMult    = 0;
        double MyTrhe = 8.0;
        for(int ipix=0;ipix<HitNPixels[HitNbCount];ipix++) {
	  if(PixelNbCount+1 > NMaxPixel) {
            cout << endl;
            cout << "Number of pixels with value " << PixelNbCount+1 << " goes beyong the maximum allowed value " << NMaxPixel << ". Exiting now!!!" << endl;
            cout << endl;
            assert(false);
          }
	
	
	  int col = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].col;
	  int row = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].row;
	  float u = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].PosUVmm.x();
	  float v = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].PosUVmm.y();

	  //Fill Branche APixelHitTAF
	  PixelHitIdx[PixelNbCount]       = HitNbCount;
	  PixelGlobalIdx[PixelNbCount]    = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].globalIdx;
	  PixelColumn[PixelNbCount]       = col;
	  PixelRow[PixelNbCount]          = row;
	  PixelAnalogCharge[PixelNbCount] = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].ChargeAnalog;
	  PixelSensorID[PixelNbCount]     = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].sensorID;
	  PixelUmm[PixelNbCount]          = u;
	  PixelVmm[PixelNbCount]          = v;
	  PixelStatus[PixelNbCount]       = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].status;

	  if(ListOfSensitivePlanes[idx_SD].digit->GetAnalysisMode() == 2 && PixelAnalogCharge[PixelNbCount] > MyTrhe*ListOfSensitivePlanes[idx_SD].digit->GetNoiseElectrons()) MyMult++;
	
	  if(row < Rowmin) Rowmin = row;
	  if(row > Rowmax) Rowmax = row;
	  if(col < Colmin) Colmin = col;
	  if(col > Colmax) Colmax = col;
	
	  double weighting_factor = 1.0;
	  if(ListOfSensitivePlanes[idx_SD].digit->GetAnalysisMode() == 2) weighting_factor = PixelAnalogCharge[PixelNbCount];
	
	  SumPosU      +=   u * weighting_factor;
	  SumPosV      +=   v * weighting_factor;
	  SumPosUU     += u*u * weighting_factor;
	  SumPosVV     += v*v * weighting_factor;
	  SumPosUV     += u*v * weighting_factor;
	
	  SumOfWeights += weighting_factor;
	
	  PixelNbCount++;
        } //End Loop on pixelMap
        SumOfWeights = TMath::Abs(SumOfWeights);
        if(ListOfSensitivePlanes[idx_SD].digit->GetAnalysisMode() == 3) MyMult = HitNPixels[HitNbCount];

        if(HitNPixels[HitNbCount] > 0 && TMath::Abs(SumOfWeights) > 1.0e-8) {
	  SumPosU  /= SumOfWeights;
	  SumPosV  /= SumOfWeights;
	  SumPosUU /= SumOfWeights;
	  SumPosUU -= pow(SumPosU,2);
	  SumPosVV /= SumOfWeights;
	  SumPosVV -= pow(SumPosV,2);
	  SumPosUV /= SumOfWeights;
	  SumPosUV -= SumPosU*SumPosV;

          HitClusterSizeCol[HitNbCount]   = Colmax - Colmin +1;
          HitClusterSizeRow[HitNbCount]   = Rowmax - Rowmin +1;
	
	  HitRecoUmm[HitNbCount]          = SumPosU;
	  HitRecoVmm[HitNbCount]          = SumPosV;
	
          HitRecoULaddermm[HitNbCount]    = -999.0;
          HitRecoVLaddermm[HitNbCount]    = -999.0;

	  float TheRMSMax;
	  float ThePhiPrimMax;
	  if(HitNPixels[HitNbCount] > 1){
	    float PhiRms = 0.5*GetATan(2*SumPosUV, SumPosUU - SumPosVV);
	  
	    float PhiRms_1;
	    float PhiRms_2;
	    if(PhiRms >= 0) {
	      PhiRms_1 = PhiRms;
	      PhiRms_2 = PhiRms - 0.5*TMath::Pi();
	    }
	    else {  
	      PhiRms_1 = PhiRms;
	      PhiRms_2 = PhiRms + 0.5*TMath::Pi();
	    }
	  
	    float sigmaU_temp1 = SumPosUU*pow(cos(PhiRms_1),2) + SumPosVV*pow(sin(PhiRms_1),2) + SumPosUV*sin(2*PhiRms_1);
	    float sigmaU_temp2 = SumPosUU*pow(cos(PhiRms_2),2) + SumPosVV*pow(sin(PhiRms_2),2) + SumPosUV*sin(2*PhiRms_2);
	    if(sigmaU_temp1 > sigmaU_temp2){
	      TheRMSMax     = sigmaU_temp1; 
	      ThePhiPrimMax = PhiRms_1;
	    }
	    else{
	      TheRMSMax     = sigmaU_temp2; 
	      ThePhiPrimMax = PhiRms_2; 
	    }
	  }
	  else{
	    TheRMSMax     = 0;
	    ThePhiPrimMax = 0;
	  }

	  HitRMSPrincipalAxis[HitNbCount] = sqrt(TMath::Abs(TheRMSMax));
	  HitPhiPrincipalAxis[HitNbCount] = ThePhiPrimMax*(180/TMath::Pi());

	  bool MyVerbose_hits = false;
	  MyVerbose_hits      = true;
	  double Limit_Delta = 30.0;
	  if((1000*TMath::Abs(HitRecoUmm[HitNbCount] - HitposAVRmmU[HitNbCount]) > Limit_Delta || 1000*TMath::Abs(HitRecoVmm[HitNbCount] - HitposAVRmmV[HitNbCount]) > Limit_Delta) && MyVerbose_hits) {
	    double Momentum  = pow(HitmomMeVX[HitNbCount],2);
	    Momentum        += pow(HitmomMeVY[HitNbCount],2);
	    Momentum        += pow(HitmomMeVZ[HitNbCount],2);
	    Momentum         = sqrt(Momentum);
	    double TheEpi    = TMath::Abs(HitposINmmW[HitNbCount] - HitposOUTmmW[HitNbCount])*1000;
	    double MyCol,MyRow;
	    double MyEdepElec = HitGeant4EdepMeV[HitNbCount]*1.0e+6/3.6;
	    ListOfSensitivePlanes[idx_SD].digit->ComputePixelPositionColRow_FromUV(HitposAVRmmU[HitNbCount]*1000,HitposAVRmmV[HitNbCount]*1000,MyCol, MyRow);
	    cout << "pos-(u,v) = (" << 1000*HitposAVRmmU[HitNbCount]                            << "," << 1000*HitposAVRmmV[HitNbCount] << ") um, "
	         << "rec-(u,v) = (" << 1000*HitRecoUmm[HitNbCount]                              << "," << 1000*HitRecoVmm[HitNbCount]   << ") um, "
	         << "del-(u,v) = (" << 1000*(HitposAVRmmU[HitNbCount] - HitRecoUmm[HitNbCount]) << "," << 1000*(HitposAVRmmV[HitNbCount] - HitRecoVmm[HitNbCount]) << ") um, "
	         << "pos-(col,row) = (" << MyCol << "," << MyRow << "),  "
	         << "SumWeights = " << SumOfWeights << ", "
	         << "Pixel-mult = " << MyMult << ",  "
	         << "RMS-princ-Axis = " << 1000*HitRMSPrincipalAxis[HitNbCount] << " um, "
	         << "Phi-princ-Axis = " << HitPhiPrincipalAxis[HitNbCount] << " deg, "
	         << "theta_loc = " << HitthetaLoc[HitNbCount] << " deg, "
	         << "phi_loc = " << HitphiLoc[HitNbCount] << " deg, "
	         << "Momentum = " << Momentum << " MeV/c, "
	         << "(W-in,W-out) = (" << 1000*HitposINmmW[HitNbCount] << "," << 1000*HitposOUTmmW[HitNbCount] << ") um, "
	         << "Delta-W = " << TheEpi << " um, "
	         << "Edep = " << MyEdepElec << " elec,  "
	         << "pdgID = " << ParticlepdgID[ipart] << ",  " 
	         << "trkID = " << ListOfParticlesInMimosa[ipart].trackID << ",  "
	         << "sensorID = " << HitsensorID[HitNbCount] << ",  "
	         << "Evt   = " << EvtNumber << ",  "
	         << endl;
	    
	    for(int ipix=0;ipix<HitNPixels[HitNbCount];ipix++) {
	      int global    = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].globalIdx;
	      int col       = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].col;
	      int row       = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].row;
	      double Charge = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix].ChargeAnalog;

	      cout << "(global,col,row,charge) = (" << global << "," << col << "," << row << "," << Charge << ")" << endl;
	    }
	  
          } // End of test if Reco hit is too far away from true hit
	
        } // End of if pixels of hit is > 0
        
      } // End of if doing digitization

      HitNbCount++; 
    } // End Loop on ListOfParticleHits
  } //End loop on Particle [ihit]

  HitNb   = HitNbCount;

  //Fill the pixels from noise
  for(int ipix=0;ipix<int(ListOfPixelsFromNoise.size());ipix++) {
    if(PixelNbCount+1 > NMaxPixel) {
      cout << endl;
      cout << "Number of pixels with value " << PixelNbCount+1 << " goes beyong the maximum allowed value " << NMaxPixel << ". Exiting now!!!" << endl;
      cout << endl;
      assert(false);
    }
    
    PixelHitIdx[PixelNbCount]       = -1;
    PixelGlobalIdx[PixelNbCount]    = ListOfPixelsFromNoise[ipix].globalIdx;
    PixelColumn[PixelNbCount]       = ListOfPixelsFromNoise[ipix].col;
    PixelRow[PixelNbCount]          = ListOfPixelsFromNoise[ipix].row;
    PixelAnalogCharge[PixelNbCount] = ListOfPixelsFromNoise[ipix].ChargeAnalog;
    PixelSensorID[PixelNbCount]     = ListOfPixelsFromNoise[ipix].sensorID;
    PixelUmm[PixelNbCount]          = ListOfPixelsFromNoise[ipix].PosUVmm.x();
    PixelVmm[PixelNbCount]          = ListOfPixelsFromNoise[ipix].PosUVmm.y();
    PixelStatus[PixelNbCount]       = ListOfPixelsFromNoise[ipix].status;
    
    PixelNbCount++;
  
  } //End Loop on pixelMap

  PixelNb = PixelNbCount;

  return;
  
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillNonSensitiveBlock(void)
{

  // This function fill the non-sensitive volumes hits block 
  
  // Filling up the particles with their hits and their turned on pixels
  NonSensitiveParticleNb = int(ListOfParticlesInNonSensitive.size());
  if(NonSensitiveParticleNb > NMaxParticle) {
    cout << endl;
    cout << "NonSensitiveParticleNb variable with value " << NonSensitiveParticleNb << " goes beyong the maximum allowed value " << NMaxParticle << ". Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  int HitNbCount   = 0;
  for(int ipart=0;ipart<int(ListOfParticlesInNonSensitive.size());ipart++) {
    ParticleBKGType[ipart] = -1;
    NonSensitiveParticlepdgID[ipart]   = ListOfParticlesInNonSensitive[ipart].pdgID;
    NonSensitiveParticleTrkVtX[ipart]  = ListOfParticlesInNonSensitive[ipart].trkVtx.x();
    NonSensitiveParticleTrkVtY[ipart]  = ListOfParticlesInNonSensitive[ipart].trkVtx.y();
    NonSensitiveParticleTrkVtZ[ipart]  = ListOfParticlesInNonSensitive[ipart].trkVtx.z();
    NonSensitiveParticleNHits[ipart]   = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits.size();
    if(NonSensitiveParticleNHits[ipart]==0) NonSensitiveParticleFirstHitIdx[ipart] = -1;
    else                                    NonSensitiveParticleFirstHitIdx[ipart] = HitNbCount;
    
    NonSensitiveParticleSensitiveIdx[ipart] = -1;
    for(int ipart2=0;ipart2<ParticleNb;ipart2++) {
      if(ParticletrkID[ipart2] == ListOfParticlesInNonSensitive[ipart].trackID) {
	NonSensitiveParticleSensitiveIdx[ipart] = ipart2;
	break;
      }
    }

    for(int ihit=0;ihit<int(ListOfParticlesInNonSensitive[ipart].ListOfParticleHits.size());ihit++) {
      if(HitNbCount+1 > kMaxHits) {
        cout << endl;
        cout << "Number of hits in non-sensitive volumes with value " << HitNbCount+1 << " goes beyong the maximum allowed value " << kMaxHits << ". Exiting now!!!" << endl;
        cout << endl;
        assert(false);
      }
      
      NonSensitiveHitParticleIdx[HitNbCount]      = ipart;
      NonSensitiveHitposINmmX[HitNbCount]         = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].posINmmXYZ.x();
      NonSensitiveHitposINmmY[HitNbCount]         = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].posINmmXYZ.y();
      NonSensitiveHitposINmmZ[HitNbCount]         = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].posINmmXYZ.z();
      NonSensitiveHitglobalTimeINns[HitNbCount]   = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].globalTime;
      NonSensitiveHitmomINMeVX[HitNbCount]        = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].momINGeV.x()*1.0e+3;
      NonSensitiveHitmomINMeVY[HitNbCount]        = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].momINGeV.y()*1.0e+3;
      NonSensitiveHitmomINMeVZ[HitNbCount]        = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].momINGeV.z()*1.0e+3;
      NonSensitiveHitEnergyINMeV[HitNbCount]      = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].EnergyINGeV*1.0e+3;
      NonSensitiveHitposOUTmmX[HitNbCount]        = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].posOUTmmXYZ.x();
      NonSensitiveHitposOUTmmY[HitNbCount]        = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].posOUTmmXYZ.y();
      NonSensitiveHitposOUTmmZ[HitNbCount]        = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].posOUTmmXYZ.z();
      NonSensitiveHitmomOUTMeVX[HitNbCount]       = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].momOUTGeV.x()*1.0e+3;
      NonSensitiveHitmomOUTMeVY[HitNbCount]       = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].momOUTGeV.y()*1.0e+3;
      NonSensitiveHitmomOUTMeVZ[HitNbCount]       = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].momOUTGeV.z()*1.0e+3;
      NonSensitiveHitEnergyOUTMeV[HitNbCount]     = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].EnergyOUTGeV*1.0e+3;
      NonSensitiveHitGeant4EdepMeV[HitNbCount]    = ListOfParticlesInNonSensitive[ipart].ListOfParticleHits[ihit].EdepMeV;
      HitNbCount++; 
    } // End Loop on ListOfParticleHits
  } //End loop on Particle [ihit]

  NonSensitiveHitNb   = HitNbCount;

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuHistoManager::FillSaturationBlock(void)
{

  // In this function will fill the saturated lines. Only works for Mi26/28 sensors.
  // Currently just letting it empty
  
  SatuNb = 0;
  //SatuLinIdx;
  //SatuSensorID;

  return;
  
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuHistoManager::FillEventSizeBlock(void)
{

  // In this function will fill be event size. It is many done for the PLUME ladders of the 
  // BEAST application (2 Ladders). Could be used to implement the full event size in one of the variables
  // i.e. EventSizeLadder1 = Actuall event size [bytes], and EventSizeLadder2 = 0
  // Currently setting it to 0
  
  EventSizeLadder1 = 0.0;
  EventSizeLadder2 = 0.0;

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::FillNtuple()
{

  // This function fills the n-tuple variables
  
  fNtuple->Fill();

  return;
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuHistoManager::EqualizeAParticleHits(AParticleHits_t &A, AParticleHits_t B)
{

  // This function equalizes two AParticleHits_t objects, A = B;
  
  A.pdgID            = B.pdgID;
  A.sensorID         = B.sensorID;
  A.ladderID         = B.ladderID;
  A.moduleID         = B.moduleID;
  A.posINmmXYZ       = B.posINmmXYZ;
  A.posINmmUVW       = B.posINmmUVW;
  A.momINGeV         = B.momINGeV;
  A.EnergyINGeV      = B.EnergyINGeV;
  A.thetaLoc         = B.thetaLoc;
  A.phiLoc           = B.phiLoc;
  A.posOUTmmXYZ      = B.posOUTmmXYZ;
  A.posOUTmmUVW      = B.posOUTmmUVW;
  A.momOUTGeV        = B.momOUTGeV;
  A.EnergyOUTGeV     = B.EnergyOUTGeV;
  A.globalTime       = B.globalTime;
  A.EdepMeV          = B.EdepMeV;
  A.clusterSizeCol   = B.clusterSizeCol;
  A.clusterSizeRow   = B.clusterSizeRow;
  A.recoUVmm         = B.recoUVmm;
  A.recoUVLaddermm   = B.recoUVLaddermm;
  A.phiPrincipalAxis = B.phiPrincipalAxis;
  A.rmsPrincipalAxis = B.rmsPrincipalAxis;
  A.thetaScattering  = B.thetaScattering;
  
  A.ListOfHitPixels.clear();
  for(int ipix=0;ipix<int(B.ListOfHitPixels.size());ipix++) A.ListOfHitPixels.push_back(B.ListOfHitPixels[ipix]);

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Long64_t MimosaSimuHistoManager::RoundOffInt(double Djentry)
{
  
  // Rounding off double to closes integer
  
  Long64_t jentry = Long64_t(Djentry);
  double delta = Djentry - jentry;
  if(delta > 0.5) jentry++;
  
  return jentry;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool  MimosaSimuHistoManager::IsAlreadyInPixelList(int Col, int Row, int TheSensorID)
{
  
  // This function checks if a pixels is already in the list of pixel turned on by particles
  
  bool IsIn = false;
  for(int ipart=0;ipart<int(ListOfParticlesInMimosa.size());ipart++) {
    for(int ihit=0;ihit<int(ListOfParticlesInMimosa[ipart].ListOfParticleHits.size());ihit++) {
      if(ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].sensorID != TheSensorID) continue;
      
      for(int ipix2=0;ipix2<int(ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels.size());ipix2++) {
	int col2 = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix2].col;
	int row2 = ListOfParticlesInMimosa[ipart].ListOfParticleHits[ihit].ListOfHitPixels[ipix2].row;
	if(Col == col2 && Row == row2) {
	  IsIn = true;
	  break;
	}
      }
      if(IsIn) break;
    }
    if(IsIn) break;
  }
    
  return IsIn;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MimosaSimuHistoManager::GetATan(double Y, double X)
{
  
  double epsilon = 1.0e-8;
  
  double delta = 0.0;
  
  if(TMath::Abs(X) < epsilon) {
    if(Y >= 0.0) delta = +0.5*TMath::Pi();
    else         delta = -0.5*TMath::Pi();
  }
  else {
    double ratio = TMath::ATan(TMath::Abs(Y/X));
    if(X > 0.0) {
      if(Y >= 0.0) delta =  ratio;
      else         delta = -ratio;
    }
    else {
      if(Y >= 0.0) delta =   TMath::Pi() - ratio;
      else         delta =  -TMath::Pi() + ratio;
    }
    
  }
  
  return delta;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuHistoManager::InitChargedParticleList(void)
{
  
  ListOfChargedParticles.clear();
  
  // Only consider stable charged particles
  
  //electron/positron
  ListOfChargedParticles.push_back(11);
  //muon+/muon-
  ListOfChargedParticles.push_back(13);
  //pi+/pi-
  ListOfChargedParticles.push_back(211);
  //K+/K-
  ListOfChargedParticles.push_back(321);
  //proton/anti-proton
  ListOfChargedParticles.push_back(2212);
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool  MimosaSimuHistoManager::IsParticleCharged(int pdgID)
{
  
  bool IsCharged = false;
  for(int ipart=0;ipart<int(ListOfChargedParticles.size());ipart++) {
    if(abs(pdgID) == ListOfChargedParticles[ipart]) {
      IsCharged = true;
      break;
    }
  }
  
  return IsCharged;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool  MimosaSimuHistoManager::IsGamma(int pdgID)
{
  
  bool IsPhoton = false;
  if(abs(pdgID) == 22) IsPhoton = true;
  
  return IsPhoton;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int  MimosaSimuHistoManager::GetIdxInSensitivePlaneList(int aSensorID)
{
  
  int idx_SD = -1;
  for(int iplane=0;iplane<int(ListOfSensitivePlanes.size());iplane++) {
    if(ListOfSensitivePlanes[iplane].sensorID == aSensorID) {
      idx_SD = iplane;
      break;
    }
  }
  
  return idx_SD;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
