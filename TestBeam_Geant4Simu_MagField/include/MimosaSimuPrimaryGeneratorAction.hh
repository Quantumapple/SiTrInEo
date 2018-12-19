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
// $Id: MimosaSimuPrimaryGeneratorAction.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MimosaSimuPrimaryGeneratorAction.hh
/// \brief Definition of the MimosaSimuPrimaryGeneratorAction class

#ifndef MimosaSimuPrimaryGeneratorAction_h
#define MimosaSimuPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "MimosaSimuSetup.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "G4RotationMatrix.hh"
#include "TH1F.h"
#include <map>

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

using namespace CLHEP;

class MimosaSimuPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MimosaSimuPrimaryGeneratorAction(MimosaSimuSetup *TheSetup);
    virtual ~MimosaSimuPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* );

    G4ParticleGun* GetParticleGun() {return fParticleGun;}
  
    // Set methods
    void   SetRandomFlag(G4bool);
    
    void   InitGenerationParameters(void);
    void   DoUniformBeamTestInit(void);
    void   DoBeamSpotBeamTestInit(void);
    void   DoSourceInit(void);
    
    void   FillParticleMap(void);
    void   CheckIfParticleIsInMap(TString AParticle);
    
    void   InitBetaSourceSpectrum(TString SourceName);
    void   InitXRaySourceSpectrum(TString SourceName);
    
    void   GenerateBeamTestPrimaries(G4Event* anEvent);
    void   GenerateSourcePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun*          fParticleGun; // G4 particle gun
    MimosaSimuSetup *fSetup;
    
    G4double nParticles;
    G4long mySeed;
    G4long myOldSeed;
    HepRandomEngine* theRandomEngine;
    RandGaussQ*      Mygauss;
    RandPoissonQ*    Mypoisson;

    TString TheParticle;
    
    //Beam-Test generation parameters
    G4double       TheBeamMomentum;
    G4ThreeVector  TheBeamDirection;
    G4double       TheBeamAngularSpreadX;
    G4double       TheBeamAngularSpreadY;
    G4double       TheBeamMomentumSpread;
    G4ThreeVector  TheBeamOrigin;
    G4double       TheBeamOriginSpreadX;
    G4double       TheBeamOriginSpreadY;
    bool           UniformOriginGeneration;
    G4ThreeVector  PerpXprime;
    G4ThreeVector  PerpYprime;
    
    //Source generation parameters
    TString                SourceType;
    G4double               SourceRadius;
    G4double               SourceHight;
    G4double               SourceActivity;
    G4double               SourceSensorROTime;
    G4ThreeVector          SourcePosition;
    G4RotationMatrix*      Source_rot;
    G4RotationMatrix*      Source_Invrot;
    G4ThreeVector          SourceLocalXVector;
    G4ThreeVector          SourceLocalYVector;
    G4ThreeVector          SourceLocalZVector;
    TH1F                   h_SourceBetaSpectrum;
    std::vector<G4double>  XRayLinesEnergy;
    std::vector<G4double>  XRayLinesProb;
    std::vector<G4double>  XRayLinesCumulProb;

    std::map<TString,G4ParticleDefinition*> particle_map;
    
    void GetPerpendicularFrame(void);
    
    bool IsBetaSource(TString SourceName);
    bool IsXRaySource(TString SourceName);
    
    G4double  Mom_beta_spectrum_pdf(G4double p,G4double Q,G4double Z,bool InclFermi);
    //G4double  Fermi_func(G4double p,G4double Q,G4double Z);
    G4double  Fermi_func(G4double p,G4double Z);
    
    G4double  GetSourceRandomMomentum(void);
    G4double  GetBetaSourceRandomMomentum(void);
    G4double  GetXRaySourceRandomMomentum(void);
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
