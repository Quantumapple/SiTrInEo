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
// $Id: MimosaSimuPrimaryGeneratorAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MimosaSimuPrimaryGeneratorAction.cc
/// \brief Implementation of the MimosaSimuPrimaryGeneratorAction class

#include "MimosaSimuPrimaryGeneratorAction.hh"

#include <assert.h>
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include <G4strstreambuf.hh>

#include "Randomize.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuPrimaryGeneratorAction::MimosaSimuPrimaryGeneratorAction(MimosaSimuSetup* TheSetup) : G4VUserPrimaryGeneratorAction(),
                                                                                                fSetup(TheSetup)
{

  mySeed       = fSetup->GetAnalysisPar().MCSeed; //time(NULL); 
  myOldSeed    = mySeed-1;
  
  theRandomEngine = new CLHEP::MTwistEngine();
  CLHEP::HepRandom::setTheEngine(theRandomEngine);
  Mygauss   = new RandGaussQ(theRandomEngine,0.0,1.0);
  Mypoisson = new RandPoissonQ(theRandomEngine,nParticles);

  theRandomEngine->setSeed(mySeed,myOldSeed);
  mySeed = G4long(theRandomEngine->flat()*100000000000.);
  
  //Filling up particle map
  FillParticleMap();
  
  fParticleGun = new G4ParticleGun(1);
  
  if(fSetup->TrackerParameter.ExpSetup == TString("Beam-Test")) {
    nParticles   = fSetup->TrackerParameter.BeamNparticles;
  
    //The user defined particle type
    TheParticle = fSetup->TrackerParameter.BeamType;
    //Check if the user defined particle type is in particle map
    CheckIfParticleIsInMap(TheParticle);
  
    TheBeamMomentum         = TMath::Abs(fSetup->TrackerParameter.BeamMomentum)*GeV;
    TheBeamDirection        = fSetup->TrackerParameter.BeamDirection;
    TheBeamAngularSpreadX   = fSetup->TrackerParameter.BeamAngularSpread(0);
    TheBeamAngularSpreadY   = fSetup->TrackerParameter.BeamAngularSpread(1);
    TheBeamMomentumSpread   = fSetup->TrackerParameter.BeamMomentumSpread;
    TheBeamOrigin           = fSetup->TrackerParameter.BeamOrigin*um;
    TheBeamOriginSpreadX    = fSetup->TrackerParameter.BeamOriginSpread(0)*um;
    TheBeamOriginSpreadY    = fSetup->TrackerParameter.BeamOriginSpread(1)*um;
  }
  else if(fSetup->CheckIfSourceExists(fSetup->TrackerParameter.ExpSetup)) {
    SourceType         = fSetup->TrackerParameter.ExpSetup;
    SourceActivity     = fSetup->TrackerParameter.SourceActivity*(1/s);
    SourceSensorROTime = fSetup->TrackerParameter.SourceSensorROTime*1.0e-6*s;
    nParticles         = 0.5*SourceActivity*SourceSensorROTime;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MimosaSimuPrimaryGeneratorAction::~MimosaSimuPrimaryGeneratorAction()
{
  if(fParticleGun) {
    delete  fParticleGun;
    fParticleGun = NULL;
  }
  if(theRandomEngine) {
    delete  theRandomEngine;
    theRandomEngine = NULL;
  }
  //delete  Mygauss;
  //delete  Mypoisson;
  
  XRayLinesEnergy.clear();
  XRayLinesProb.clear();
  XRayLinesCumulProb.clear();
  particle_map.clear();
  
  //cout << "MimosaSimuPrimaryGeneratorAction::~MimosaSimuPrimaryGeneratorAction:: All OK!!!" << endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
  
  if(anEvent->GetEventID()+1 == 1) InitGenerationParameters();

  if(fSetup->TrackerParameter.ExpSetup == TString("Beam-Test"))            GenerateBeamTestPrimaries(anEvent);
  else if(fSetup->CheckIfSourceExists(fSetup->TrackerParameter.ExpSetup))  GenerateSourcePrimaries(anEvent);

  return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   MimosaSimuPrimaryGeneratorAction::GenerateBeamTestPrimaries(G4Event* anEvent)
{
  
  bool Myverbose = false;
  //Myverbose = true;
  
  G4int TheNparticles;
  if(fSetup->TrackerParameter.BeamRandNparticles) {
    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    TheNparticles = Mypoisson->shoot(nParticles);
  }
  else TheNparticles = nParticles;
  
  if(Myverbose) { 
    G4cout << G4endl;
    G4cout << "Generating " << TheNparticles << " particles for event " << anEvent->GetEventID()+1 << " (" << nParticles << ")" << G4endl;
  }
  
  for(G4int ipar  = 0; ipar < TheNparticles; ipar++) {    
    G4ThreeVector GenerationPostion(0.0,0.0,0.0);
    G4ThreeVector MomentumDirection(0.0,0.0,1.0);
    G4double      ParticleMomentum = 0.005*GeV;
    
    mySeed++;
    if(mySeed==myOldSeed) std::cout<<"CAUTION SAME SEED AS BEFORE !!!!!!!!!!!!!!!!!!!"<<std::endl;
    myOldSeed = mySeed;
    
    //Beam-Test generation
    if(UniformOriginGeneration) {
      mySeed++;
      theRandomEngine->setSeed(mySeed,myOldSeed);
      GenerationPostion(0) = (2.0*theRandomEngine->flat() - 1.0)*TheBeamOriginSpreadX;
      
      mySeed++;
      theRandomEngine->setSeed(mySeed,myOldSeed);
      GenerationPostion(1) = (2.0*theRandomEngine->flat() - 1.0)*TheBeamOriginSpreadY;
      
      GenerationPostion(2) = TheBeamOrigin(2);
    }
    else {
      mySeed++;
      theRandomEngine->setSeed(mySeed,myOldSeed);
      G4double DeltaXprime = Mygauss->shoot()*TheBeamOriginSpreadX;
      
      mySeed++;
      theRandomEngine->setSeed(mySeed,myOldSeed);
      G4double DeltaYprime = Mygauss->shoot()*TheBeamOriginSpreadY;
      
      GenerationPostion  = TheBeamOrigin;
      GenerationPostion += DeltaXprime*PerpXprime;
      GenerationPostion += DeltaYprime*PerpYprime;
    }
      
    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    G4double fNumXAngle = Mygauss->shoot()*TheBeamAngularSpreadX;

    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    G4double fNumYAngle = Mygauss->shoot()*TheBeamAngularSpreadY;

    G4double Ppgen = 1.0/sqrt(1.0 + pow(fNumXAngle,2) + pow(fNumYAngle,2));
    G4double Pxgen = fNumXAngle*Ppgen;
    G4double Pygen = fNumYAngle*Ppgen;
      
    MomentumDirection  = Pxgen*PerpXprime;
    MomentumDirection += Pygen*PerpYprime;
    MomentumDirection += Ppgen*TheBeamDirection;
    MomentumDirection *= (1.0/MomentumDirection.mag());
      
    mySeed++;
    myOldSeed++;
  
    ParticleMomentum = (Mygauss->shoot()*TheBeamMomentumSpread + 1.0)*TheBeamMomentum;
    
    if(Myverbose) { 
      G4cout << "Generating particle " << ipar+1 << ": " << G4endl;
      G4cout << "   name      = "  << particle_map[TheParticle]->GetParticleName() << G4endl;
      G4cout << "   Position  = (" << GenerationPostion(0)/mm << "," << GenerationPostion(1)/mm << "," << GenerationPostion(2)/mm << ") mm" << G4endl;
      G4cout << "   Direction = (" << MomentumDirection(0)/mm << "," << MomentumDirection(1)/mm << "," << MomentumDirection(2)/mm << ") mm" << G4endl;
      G4cout << "   Momentum  =  " << ParticleMomentum/GeV << " GeV" << G4endl;
    }
    
    //G4ParticleMomentum TheFinalMomentum = ParticleMomentum*MomentumDirection;
    
    G4strstreambuf* oldBuffer = dynamic_cast<G4strstreambuf*>(G4cout.rdbuf(0));
    
    //fParticleGun->SetParticleDefinition(particle_map[TheParticle]);
    fParticleGun->SetParticlePosition(GenerationPostion);
    fParticleGun->SetParticleMomentumDirection(MomentumDirection);
    fParticleGun->SetParticleMomentum(ParticleMomentum);
    //fParticleGun->SetParticleMomentum(TheFinalMomentum);
    
    G4cout.rdbuf(oldBuffer);
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
  }
  if(Myverbose) G4cout << G4endl;
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   MimosaSimuPrimaryGeneratorAction::GenerateSourcePrimaries(G4Event* anEvent)
{
  
  bool Myverbose = false;
  //Myverbose = true;

  G4int TheNparticles; 
  mySeed++;
  theRandomEngine->setSeed(mySeed,myOldSeed);
  TheNparticles = Mypoisson->shoot(nParticles);

  if(Myverbose) { 
    G4cout << G4endl;
    G4cout << "Generating " << TheNparticles << " particles for event " << anEvent->GetEventID()+1 << " (" << nParticles*0.5 << ")" << G4endl;
  }
  
  for(G4int ipar  = 0; ipar < TheNparticles; ipar++) {    
    G4ThreeVector GenerationPostion(0.0,0.0,0.0);
    G4ThreeVector MomentumDirection(0.0,0.0,1.0);
    G4double      ParticleMomentum = 0.005*GeV;
    
    mySeed++;
    if(mySeed==myOldSeed) std::cout<<"CAUTION SAME SEED AS BEFORE !!!!!!!!!!!!!!!!!!!"<<std::endl;
    myOldSeed = mySeed;
    
    //Uniform distribution in source disk
    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    G4double phi_pos = theRandomEngine->flat()*2.0*TMath::Pi();
    
    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    G4double radius  = theRandomEngine->flat()*SourceRadius;
    
    GenerationPostion = G4ThreeVector(radius*cos(phi_pos),
				      radius*sin(phi_pos),
				      -0.5*SourceHight - 0.01*mm);
    GenerationPostion.transform(*Source_rot);
    //GenerationPostion.transform(*Source_Invrot);
    GenerationPostion += SourcePosition;
    
    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    G4double phi_dir = theRandomEngine->flat()*2.0*TMath::Pi();
    
    mySeed++;
    theRandomEngine->setSeed(mySeed,myOldSeed);
    G4double costheta_dir = theRandomEngine->flat();
    
    MomentumDirection  = cos(phi_dir)*sqrt(1.0 - pow(costheta_dir,2))*SourceLocalXVector;
    MomentumDirection += sin(phi_dir)*sqrt(1.0 - pow(costheta_dir,2))*SourceLocalYVector;
    MomentumDirection +=                                 costheta_dir*SourceLocalZVector;
    
    ParticleMomentum = GetSourceRandomMomentum();
    
    if(Myverbose) {
      G4cout << "Generating particle " << ipar+1 << ": " << G4endl;
      G4cout << "   name      = "  << particle_map[TheParticle]->GetParticleName() << G4endl;
      G4cout << "   Position  = (" << GenerationPostion(0)/mm << "," << GenerationPostion(1)/mm << "," << GenerationPostion(2)/mm << ") mm" << G4endl;
      G4cout << "   Direction = (" << MomentumDirection(0)/mm << "," << MomentumDirection(1)/mm << "," << MomentumDirection(2)/mm << ") mm" << G4endl;
      G4cout << "   Momentum  =  " << ParticleMomentum/MeV << " MeV" << G4endl;
    }
    
    G4strstreambuf* oldBuffer = dynamic_cast<G4strstreambuf*>(G4cout.rdbuf(0));
    
    //fParticleGun->SetParticleDefinition(particle_map[TheParticle]);
    fParticleGun->SetParticlePosition(GenerationPostion);
    fParticleGun->SetParticleMomentumDirection(MomentumDirection);
    fParticleGun->SetParticleMomentum(ParticleMomentum);
    
    G4cout.rdbuf(oldBuffer);
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
  }
  if(Myverbose) G4cout << G4endl;
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::InitGenerationParameters(void)
{

  if(fSetup->TrackerParameter.ExpSetup == TString("Beam-Test")) {
    //Beam-Test experimental setup
    UniformOriginGeneration = false;
    if(fSetup->TrackerParameter.BeamOrigin == G4ThreeVector(0.0,0.0,-1.0e+6) || fSetup->TrackerParameter.BeamOriginSpread(0) < 0.0 || fSetup->TrackerParameter.BeamOriginSpread(1) < 0.0) UniformOriginGeneration = true;
    
    if(UniformOriginGeneration) DoUniformBeamTestInit();
    else                        DoBeamSpotBeamTestInit();
    
    G4cout << "Mean energy of " << TheBeamMomentum/GeV << " GeV, with smearing of " << TheBeamMomentumSpread*100.0 << " %" << G4endl;
    G4cout << "Mean angular dispersion RMS in (X',Y') = (" << TheBeamAngularSpreadX << "," << TheBeamAngularSpreadY << ")" << G4endl;
    G4cout << G4endl;
  }
  else if(fSetup->CheckIfSourceExists(fSetup->TrackerParameter.ExpSetup)) {
    //Source initiaization
    DoSourceInit(); 
  }
  
  fParticleGun->SetParticleDefinition(particle_map[TheParticle]);

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::DoUniformBeamTestInit(void)
{

  G4cerr << "BeamOrigin or BeamOriginSpread  parameters not specified. Will perform uniform generation in lowest Z face of the world volume cube" << G4endl;
    
  G4double worldZHalfLength = 0;
  G4double worldXHalfLength = 0;
  G4double worldYHalfLength = 0;
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = NULL;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();
    worldXHalfLength = worldBox->GetXHalfLength();
    worldYHalfLength = worldBox->GetYHalfLength();
  }
  else  {
    G4cerr << "World volume of box not found."       << G4endl;
    G4cerr << "Perhaps you have changed geometry."   << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
    assert(false);
  }

  // Z-vertex position: 
  TheBeamOrigin(0) = 0.0;
  TheBeamOrigin(0) = 1.0;
  TheBeamOrigin(2) = -worldZHalfLength + 0.001*worldZHalfLength;
  if(TheBeamDirection.dot(G4ThreeVector(0.0,0.0,1.0)) < 0.0) TheBeamOrigin(2) = +worldZHalfLength - 0.001*worldZHalfLength;

  //Range in X: (Flat distribution)
  TheBeamOriginSpreadX = worldXHalfLength*(1.0 - 0.001);
  TheBeamOriginSpreadY = worldYHalfLength*(1.0 - 0.001);
  
  G4cout << G4endl;
  G4cout << "Beam particles will be generated with " << G4endl;
  G4cout << "Z-origin = " << TheBeamOrigin(2)/mm << " mm,"  << G4endl;
  G4cout << "X and Y origin with flat distribution in XRange x YRange of ("
         << -TheBeamOriginSpreadX/mm << "," << TheBeamOriginSpreadX/mm << ") x (" << -TheBeamOriginSpreadX/mm << "," << TheBeamOriginSpreadX/mm << ") mm^2"
         << G4endl;
  G4cout << G4endl;

  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::DoBeamSpotBeamTestInit(void)
{

  GetPerpendicularFrame();
    
  G4cout << G4endl;
  G4cout << "The beam direction = (" << TheBeamDirection(0) << "," << TheBeamDirection(1) << "," << TheBeamDirection(2) << ")" << G4endl;
  G4cout << "X' director vector = (" << PerpXprime(0)       << "," << PerpXprime(1)       << "," << PerpXprime(2)       << ")" << G4endl;
  G4cout << "Y' director vector = (" << PerpYprime(0)       << "," << PerpYprime(1)       << "," << PerpYprime(2)       << ")" << G4endl;
  G4cout << " beam direction dot X' dir = " << TheBeamDirection.dot(PerpXprime) << endl;
  G4cout << " beam direction dot Y' dir = " << TheBeamDirection.dot(PerpYprime) << endl;
  G4cout << " X' dir dot Y' dir         = " << PerpXprime.dot(PerpYprime) << endl;
  G4cout << G4endl;
  
  G4cout << G4endl;
  G4cout << "Beam particles will be generated at origin = ( " << TheBeamOrigin(0)/mm << "," << TheBeamOrigin(1)/mm << "," << TheBeamOrigin(2)/mm << ") mm" << G4endl;
  G4cout << " with (X',Y') spread of ( " << TheBeamOriginSpreadX/mm << "," << TheBeamOriginSpreadY/mm << ") mm" << endl;
  G4cout << G4endl;

  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::DoSourceInit(void)
{
  
  //Getting the source Logical volume
  G4LogicalVolume* sourceLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Source_LV");
  G4Tubs* sourceTube = NULL;
  if ( sourceLV ) sourceTube = dynamic_cast<G4Tubs*>(sourceLV->GetSolid());
  if ( sourceTube ) {
    //Get the radius and hight
    SourceRadius =     sourceTube->GetOuterRadius();
    SourceHight  = 2.0*sourceTube->GetZHalfLength();
  }
  else  {
    G4cerr << "Source_LV volume of Tubs not found."  << G4endl;
    G4cerr << "Perhaps you have changed geometry."   << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
    assert(false);
  }
  
  //Getting the source Physical volume
  G4VPhysicalVolume* sourcePV = G4PhysicalVolumeStore::GetInstance()->GetVolume("Source_PV");
  if (sourcePV) {
    //Get the volume position
    SourcePosition = sourcePV->GetObjectTranslation();
    
    //Get the volume rotation matrix
    Source_rot     = new G4RotationMatrix();
    Source_Invrot  = new G4RotationMatrix();
    *Source_rot    = sourcePV->GetObjectRotationValue();
    *Source_Invrot = Source_rot->inverse();
  }
  else  {
    G4cerr << "Source_PV volume not found."          << G4endl;
    G4cerr << "Perhaps you have changed geometry."   << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
    assert(false);
  }
  
  G4RotationMatrix* tmp_rot = Source_rot;
  //G4RotationMatrix* tmp_rot = Source_rot;
  
  SourceLocalXVector = G4ThreeVector(1.0,0.0,0.0);
  SourceLocalXVector.transform(*tmp_rot);
  SourceLocalXVector *= -1.0;
  
  SourceLocalYVector = G4ThreeVector(0.0,1.0,0.0);
  SourceLocalYVector.transform(*tmp_rot);
  SourceLocalYVector *= -1.0;
  
  SourceLocalZVector = G4ThreeVector(0.0,0.0,1.0);
  SourceLocalZVector.transform(*tmp_rot);
  SourceLocalZVector *= -1.0;
  
  //Initializing momentum or energy spectrum of sources
  if(IsBetaSource(SourceType))       InitBetaSourceSpectrum(SourceType);
  else if(IsXRaySource(SourceType))  InitXRaySourceSpectrum(SourceType);
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::GetPerpendicularFrame(void)
{

  G4double AbsPx = TMath::Abs(TheBeamDirection(0));
  G4double AbsPy = TMath::Abs(TheBeamDirection(1));
  G4double AbsPz = TMath::Abs(TheBeamDirection(2));
    
  G4double MinPcompos = TMath::Max(AbsPx,AbsPy);
  MinPcompos = TMath::Max(AbsPz,MinPcompos);
    
  if(TMath::Abs(MinPcompos - AbsPz) < 1.0e-8) {
    //p ~ z
    //X' = X & Y' = Y
    G4cout << " p ~ z, X' = X & Y' = Y" << endl;
    
    G4double b = 0.0;
    G4double a = 1.0/sqrt(1.0 + pow(AbsPx/AbsPz,2));
    G4double c = -(AbsPx/AbsPz)*a;
    
    PerpXprime = G4ThreeVector(a,b,c);
  }
  else if(TMath::Abs(MinPcompos - AbsPx) < 1.0e-8) {
    //p ~ x
    //X' = Z & Y' = Y
    G4cout << " p ~ x, X' = Z & Y' = Y" << endl;
    
    G4double b = 0.0;
    G4double c = 1.0/sqrt(1.0 + pow(AbsPz/AbsPx,2));
    G4double a = -(AbsPz/AbsPx)*c;
    
    PerpXprime = G4ThreeVector(a,b,c);
  }
  else if(TMath::Abs(MinPcompos - AbsPy) < 1.0e-8) {
    //p ~ y
    //X' = X & Y' = Z
    G4cout << " p ~ y, X' = X & Y' = Z" << endl;
    
    G4double b = 0.0;
    G4double a = 1.0/sqrt(1.0 + pow(AbsPx/AbsPy,2));
    G4double c = -(AbsPx/AbsPy)*a;
    
    PerpXprime = G4ThreeVector(a,b,c);
  }
  
  PerpYprime = TheBeamDirection.cross(PerpXprime);
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuPrimaryGeneratorAction::FillParticleMap(void)
{
  
  // default particle kinematic
  G4ParticleTable* particleTable            = G4ParticleTable::GetParticleTable();
  //The particle map
  particle_map[TString("electron")]    = particleTable->FindParticle("e-");
  particle_map[TString("positron")]    = particleTable->FindParticle("e+");
  particle_map[TString("pion-")]       = particleTable->FindParticle("pi-");
  particle_map[TString("pion+")]       = particleTable->FindParticle("pi+");
  particle_map[TString("Kaon-")]       = particleTable->FindParticle("kaon-");
  particle_map[TString("Kaon+")]       = particleTable->FindParticle("kaon+");
  particle_map[TString("proton")]      = particleTable->FindParticle("proton");
  particle_map[TString("anti-proton")] = particleTable->FindParticle("anti_proton");
  particle_map[TString("gamma")]       = particleTable->FindParticle("gamma");
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuPrimaryGeneratorAction::CheckIfParticleIsInMap(TString AParticle)
{
  
  bool IsInMap = false;
  std::map<TString, G4ParticleDefinition*>::iterator iter;
  for(iter = particle_map.begin(); iter != particle_map.end(); iter++) {
    TString ParNameTmp = iter->first;
    if(AParticle == ParNameTmp) {
      IsInMap = true;
      break;
    }
  }
  if(!IsInMap) {
    cout << endl;
    cout << "The specified particle type " << TheParticle.Data() << " is not in internal database:" << endl;
    std::map<TString, G4ParticleDefinition*>::iterator iter2;
    for(iter2 = particle_map.begin(); iter2 != particle_map.end(); iter2++) {
      TString ParNameTmp = iter2->first;
      cout << " - " << ParNameTmp.Data() << endl;
    }
    cout << "Check your inputs. Exciting now!!!" << endl;
    cout << endl;
    assert(false);
  }
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool MimosaSimuPrimaryGeneratorAction::IsBetaSource(TString SourceName)
{
  
  bool Isbeta_source = false;
  if(SourceName == TString("Source-Sr90")) Isbeta_source = true;
  
  return Isbeta_source;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool MimosaSimuPrimaryGeneratorAction::IsXRaySource(TString SourceName)
{
  
  bool Isxray_source = false;
  if(SourceName == TString("Source-Fe55")) Isxray_source = true;
  
  return Isxray_source;
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuPrimaryGeneratorAction::InitBetaSourceSpectrum(TString SourceName)
{
  
  if(SourceName == TString("Source-Sr90")) {
    TheParticle = TString("electron");
    int Z_dau_Sr_90 = 39;
    int Z_dau_Y_90  = 40;

    double Q_Sr_90       = 0.546*MeV;
    double p_Sr_90_limit = sqrt(pow(Q_Sr_90,2) - pow(particle_map[TheParticle]->GetPDGMass(),2));
    double Q_Y_90        = 2.280*MeV;
    double p_Y_90_limit  = sqrt(pow(Q_Y_90,2) - pow(particle_map[TheParticle]->GetPDGMass(),2));
    
    double porcent = 0.10;
    porcent = 0.0;
    int Nbins      = 2000;
    
    double R_mom_Y_90[2];
    R_mom_Y_90[0] = 0.0;
    R_mom_Y_90[1] = p_Y_90_limit + porcent*(p_Y_90_limit - 0.0);
    
    double R_mom_Sr_90[2];
    R_mom_Sr_90[0] = 0.0;
    R_mom_Sr_90[1] = p_Sr_90_limit + porcent*(p_Sr_90_limit - 0.0);
    
    double R_mom_SrY_90[2];
    R_mom_SrY_90[0] = 0.0;
    R_mom_SrY_90[1] = TMath::Max(R_mom_Sr_90[1],R_mom_Y_90[1]);
    
    //cout << "R_mom_Sr_90  (" << R_mom_Sr_90[0]/MeV  << "," << R_mom_Sr_90[1]/MeV  << ") MeV" << endl;
    //cout << "R_mom_Y_90   (" << R_mom_Y_90[0]/MeV   << "," << R_mom_Y_90[1]/MeV   << ") MeV" << endl;
    //cout << "R_mom_SrY_90 (" << R_mom_SrY_90[0]/MeV << "," << R_mom_SrY_90[1]/MeV << ") MeV" << endl;

    TH1F h_MomSpec_Sr_90("h_MomSpec_Sr_90",
			 "#beta^{-} momentum spectrum from ^{90}_{38}Sr #rightarrow ^{90}_{39}Y + e^{-} + #bar{#nu}_{e}",
			 Nbins,R_mom_SrY_90[0],R_mom_SrY_90[1]);
    TH1F h_MomSpec_Y_90("h_MomSpec_Y_90",
			"#beta^{-} momentum spectrum from ^{90}_{39}Sr #rightarrow ^{90}_{39}Y #rightarrow ^{90}_{40}Zr + e^{-} + #bar{#nu}_{e}",
			Nbins,R_mom_SrY_90[0],R_mom_SrY_90[1]);

    h_SourceBetaSpectrum = TH1F("h_SourceBetaSpectrum",
				"p spectrum in T_{0} from ^{90}_{38}Sr and ^{90}_{39}Y #beta decays",
				Nbins,R_mom_SrY_90[0],R_mom_SrY_90[1]);

    for(int i=0;i<Nbins;i++) {
      double val,p_val;
      
      p_val = R_mom_SrY_90[0] + (i+0.5)*(R_mom_SrY_90[1] - R_mom_SrY_90[0])/Nbins;
      
      val   = Mom_beta_spectrum_pdf(p_val,Q_Sr_90,Z_dau_Sr_90,true);
      h_MomSpec_Sr_90.SetBinContent(i+1,val);
      
      val = Mom_beta_spectrum_pdf(p_val,Q_Y_90,Z_dau_Y_90,true);
      h_MomSpec_Y_90.SetBinContent(i+1,val);

    }
    h_MomSpec_Sr_90.Scale(1.0/h_MomSpec_Sr_90.Integral("width"));
    h_MomSpec_Y_90.Scale(1.0/h_MomSpec_Y_90.Integral("width"));

    for(int i=0;i<Nbins;i++) {
      double val;
      //double p_val = R_mom_SrY_90[0] + (i+0.5)*(R_mom_SrY_90[1] - R_mom_SrY_90[0])/Nbins;
      
      val  = 0.5*h_MomSpec_Sr_90.GetBinContent(i+1);
      val += 0.5*h_MomSpec_Y_90.GetBinContent(i+1);
      h_SourceBetaSpectrum.SetBinContent(i+1,val);
    }

  }
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuPrimaryGeneratorAction::InitXRaySourceSpectrum(TString SourceName)
{
  
  TheParticle = TString("gamma");
  
  XRayLinesEnergy.clear();
  XRayLinesProb.clear();
  XRayLinesCumulProb.clear();
  
  if(SourceName == TString("Source-Fe55")) {
    G4double Energy,Prob;
    
    Energy = 5.9*keV;
    Prob   = 0.882;
    XRayLinesEnergy.push_back(Energy);
    XRayLinesProb.push_back(Prob);
    
    Energy = 6.5*keV;
    Prob   = 0.118;
    XRayLinesEnergy.push_back(Energy);
    XRayLinesProb.push_back(Prob);
  }
  
  G4double CumulProb = 0.0;
  for(int i=0;i<int(XRayLinesProb.size());i++) {
    CumulProb += XRayLinesProb[i];
    XRayLinesCumulProb.push_back(CumulProb);
  }
  
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MimosaSimuPrimaryGeneratorAction::Mom_beta_spectrum_pdf(G4double p,
							         G4double Q,
							         G4double Z,
							         bool InclFermi)
{
  
  G4double E = sqrt(pow(p,2) + pow(particle_map[TheParticle]->GetPDGMass(),2));
  if(E > Q)   return 0.0;
  if(p < 0.0) return 0.0;

  G4double pdf   = pow(p,2)*pow(Q-E,2);
  //if(InclFermi) pdf *= Fermi_func(p,Q,Z);
  if(InclFermi) pdf *= Fermi_func(p,Z);

  return pdf;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MimosaSimuPrimaryGeneratorAction::Fermi_func(G4double p,
						      //G4double Q,
						      G4double Z)
{

  G4double alpha = 7.29735257e-3;
  
  G4double E = sqrt(pow(p,2) + pow(particle_map[TheParticle]->GetPDGMass(),2));
  //G4double S   = sqrt(1.0 - pow(alpha*Z,2));
  //G4double eta = -Z*alpha*E/p;

  G4double a = 2*TMath::Pi()*alpha*Z;
  G4double b = a/(1.0 - TMath::Exp(-a));
  G4double C = b - a;
  G4double d = 0.5*(b - 1.0);

  G4double val = a*(E/p) + C/(1.0 + d/pow(p,2));

  return val;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double  MimosaSimuPrimaryGeneratorAction::GetSourceRandomMomentum(void)
{
  
  G4double Momentum = 0.0;
  
  if(IsBetaSource(SourceType))      Momentum = GetBetaSourceRandomMomentum();
  else if(IsXRaySource(SourceType)) Momentum = GetXRaySourceRandomMomentum();
  else Momentum = 0.0;
  
  return Momentum;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double  MimosaSimuPrimaryGeneratorAction::GetBetaSourceRandomMomentum(void)
{
  
  G4double Momentum = 0.0;
  
  Momentum = h_SourceBetaSpectrum.GetRandom();
 
  return Momentum;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MimosaSimuPrimaryGeneratorAction::GetXRaySourceRandomMomentum()
{
  
  G4double Momentum = 1.0e-6;
  
  mySeed++;
  theRandomEngine->setSeed(mySeed,myOldSeed);
  G4double flat_num = theRandomEngine->flat();
  
  if(XRayLinesCumulProb.size() == 1) Momentum = XRayLinesEnergy[0];
  else {
    for(int i=0;i<int(XRayLinesCumulProb.size());i++) {
      G4double Prob1,Prob2;
      if(i==0) Prob1 = 0.0;
      else     Prob1 = XRayLinesCumulProb[i-1];
      Prob2 = XRayLinesCumulProb[i];
      if(flat_num >= Prob1 && flat_num < Prob2) {
        Momentum = XRayLinesEnergy[i];
        break;
      }
    }
  }
  
  return Momentum;
  
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

