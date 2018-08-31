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
/// \file analysis/AnaEx02/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
// $Id: HistoManager.hh 74272 2013-10-02 14:48:50Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MimosaSimuHistoManager_h
#define MimosaSimuHistoManager_h 1

#include "MimosaSimuTrackerHit.hh"
#include "MimosaSimuTrajectoryHit.hh"
#include "MimosaSimuSetup.hh"

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TRandom.h"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include "digit_b2.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 class TFile;
 class TTree;
 class TH1D;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int kMaxHits(4000);

const G4int NMaxParticle(300);  //Maximum number of particle       save in the fNTuple
const G4int NMaxHit(3000);      //Maximum number of hit            save in the fNTuple
const G4int NMaxPixel(50000000);  //Maximum number of pixel          save in the fNTuple
const G4int NMaxSatu(1000);     //Maximum number of Saturated line save in the fNTuple

class MimosaSimuHistoManager
{
  
  private:

    // Geant4 hit structure
    struct Hit_t {
      G4int            pdgID;         // particle pdg ID
      G4int            sensorID;      // sensorID
      G4int            trackID;       // particle track ID
      G4ThreeVector    trackVtx;      // particle origin vertex 
      G4double         energyDepMeV;  // deposited energy by ionization, in MeV
      G4ThreeVector    posINmm_XYZ;   // enter position in volume in global ref. frame (XYZ), in mm
      G4ThreeVector    posINmm_UVW;   // enter position in volume in local  ref. frame (UVW), in mm
      G4double         thetaLoc;      // enter theta angle, in deg
      G4double         phiLoc;        // enter phi   angle, in deg
      G4LorentzVector  momINGeV;      // enter 3-momentum, in GeV
      G4double         EnergyINGeV;   // enter Energy, in GeV
      G4double         globalTimeIN;  // enter global time, in ns
      G4ThreeVector    posOUTmm_XYZ;  // exit  position in volume in global ref. frame (XYZ), in mm
      G4ThreeVector    posOUTmm_UVW;  // exit  position in volume in local  ref. frame (UVW), in mm
      G4LorentzVector  momOUTGeV;     // exit  3-momentum, in GeV
      G4double         EnergyOUTGeV;  // exit  Energy, in GeV
      G4double         globalTimeOUT; // exit  global time, in ns
    };
    Hit_t MyHit;
    std::vector<Hit_t>  HitsListInMimosa;        // List of hits in mimosa (sensitive volumes)
    std::vector<Hit_t>  HitsListInNonSensitive;  // List of hits in non-sensitive volumes
    
    // structure of hit pixels (fill after digitization)
    struct AHitPixels_t {
      G4int         globalIdx;    // pixel global index
      G4int         col;          // pixel column
      G4int         row;          // pixel row
      G4double      ChargeAnalog; // pixel collected charge, in electrons (ADC) for digital (analogue) output sensis, i.e. AnalysisMode = 3 (2)
      G4int         sensorID;     // sensor ID
      G4ThreeVector PosUVmm;      // pixel position (U,V), in mm
      G4int         status;       // status variabe. It is always 1 at least the pixel info is lost by saturation, in which case is 0
    };
    std::vector<AHitPixels_t> ListOfPixelsFromNoise;

    // structure of hits produced by a particle
    struct AParticleHits_t {
      G4int          pdgID;                       // pdg ID of particle producing hit
      G4int          sensorID;                    // sensor ID
      G4int          ladderID;                    // ladder ID, set always to -1. variable used for PLUME application BEAST
      G4int          moduleID;                    // module ID, set always to -1. variable used for PLUME application BEAST
      G4ThreeVector  posINmmXYZ;                  // enter position in volume in global ref. frame (XYZ), in mm
      G4ThreeVector  posINmmUVW;                  // enter position in volume in local  ref. frame (XYZ), in mm
      G4ThreeVector  momINGeV;                    // enter 3-momentum, in GeV
      G4double       EnergyINGeV;                 // enter Energy, in GeV
      G4double       thetaLoc;                    // enter theta angle, in deg
      G4double       phiLoc;                      // enter theta phi,   in deg
      G4ThreeVector  posOUTmmXYZ;                 // exit  position in volume in global ref. frame (XYZ), in mm
      G4ThreeVector  posOUTmmUVW;                 // exit  position in volume in local  ref. frame (XYZ), in mm
      G4ThreeVector  momOUTGeV;                   // enter 3-momentum, in GeV
      G4double       EnergyOUTGeV;                // enter Energy, in GeV
      G4double       globalTime;                  // enter global time, in ns
      G4double       EdepMeV;                     // deposited energy by ionization, in MeV
      G4double       clusterSizeCol;              // cluster size in column direction => only set if digitization turned on
      G4double       clusterSizeRow;              // cluster size in row    direction => only set if digitization turned on
      G4ThreeVector  recoUVmm;                    // reconstructed local hit position (U,V), in mm => only set if digitization turned on
      G4ThreeVector  recoUVLaddermm;              // variable for PLUME BEAST application, always set ot -999.0
      G4double       phiPrincipalAxis;            // phi angle of cluster principal axis, in deg
      G4double       rmsPrincipalAxis;            // rms of cluster principal axis, in um
      G4double       thetaScattering;             // always set to -1
      std::vector<AHitPixels_t> ListOfHitPixels;  // list of hit pixels => only filled up if digitization turned on
    };

    // structure with particle variables
    struct AParticle_t {
      G4int          pdgID;                             // particle pdg ID
      G4int          trackID;                           // particle track ID
      G4ThreeVector  trkVtx;                            // particle origin vertex
      G4double       globalTime;                        // particle origin global time
      std::vector<AParticleHits_t> ListOfParticleHits;  // list of hits produced by this particle
    };
    std::vector<AParticle_t>  ListOfParticlesInMimosa;       // List of particles producing hits in mimosa (sensitive volumes)
    std::vector<AParticle_t>  ListOfParticlesInNonSensitive; // List of particles producing hits non-sensitive volumes
    
    // structed with variables of a sensitive plane (Readout > 0)
    // variables contain all the information for the digitization
    struct ASensitivePlane_t {
      G4int      sensorID;       // sensor ID
      float      DigitizeOcc;    // Occupancy level used to fill pixel from noise list
      float      DigitizeThre;   // Threshold level (in noise units), to perform digitization only for digital output sensors 
      
      DIGIT_B2*  digit;          // digitizer object. It will perform the digitization
    };
    std::vector<ASensitivePlane_t> ListOfSensitivePlanes; // List with all the sensitive planes
    
    TRandom *rand;   // random number generator object

    std::vector<int>  ListOfChargedParticles;  //List of pdg ID of stable charged particles
    
    MimosaSimuSetup* fSetup; // pointer to Setup object
    TFile*   fRootFile;
    TTree*   fNtuple;
    
    //Mimosa Hit variables
    //Definition of Particles Branch
    Int_t    ParticleNb;                          // number of particles hitting mimosa (sensitive volume) in event
    Int_t    ParticleBKGType[NMaxParticle];       // varibales always set -1. mainly used for PLUME BEAST application
    Int_t    ParticletrkID[NMaxParticle];         // particle trk ID
    Int_t    ParticlepdgID[NMaxParticle];         // particle pdg ID
    Float_t  ParticleTrkVtX[NMaxParticle];        // particle's origin vertex, X coordinate, in mm
    Float_t  ParticleTrkVtY[NMaxParticle];        // particle's origin vertex, Y coordinate, in mm
    Float_t  ParticleTrkVtZ[NMaxParticle];        // particle's origin vertex, Z coordinate, in mm
    Int_t    ParticleNHits[NMaxParticle];         // particle's number of hits
    Int_t    ParticleFirstHitIdx[NMaxParticle];   // index on hit block of particle's 1st hit
   
    //Definition of Hits Branch
    Int_t   HitNb;                         // Number of hits in mimosa in event
    Int_t   HitParticleIdx[NMaxHit];       // index in particle block of particle producing this hit
   
    Int_t   HitsensorID[NMaxHit];          // sensir ID
    Int_t   HitladderID[NMaxHit];          // always set to -1
    Int_t   HitmoduleID[NMaxHit];          // always set to -1
   
    Float_t HitposINmmX[NMaxHit];          // enter position in global ref. frame, X coordinate, in mm
    Float_t HitposINmmY[NMaxHit];          // enter position in global ref. frame, Y coordinate, in mm 
    Float_t HitposINmmZ[NMaxHit];          // enter position in global ref. frame, Z coordinate, in mm
   
    Float_t HitposINmmU[NMaxHit];          // enter position in local  ref. frame, U coordinate, in mm
    Float_t HitposINmmV[NMaxHit];          // enter position in local  ref. frame, V coordinate, in mm
    Float_t HitposINmmW[NMaxHit];          // enter position in local  ref. frame, W coordinate, in mm
  
    Float_t HitposOUTmmU[NMaxHit];         // exit  position in local  ref. frame, U coordinate, in mm
    Float_t HitposOUTmmV[NMaxHit];         // exit  position in local  ref. frame, V coordinate, in mm
    Float_t HitposOUTmmW[NMaxHit];         // exit  position in local  ref. frame, W coordinate, in mm
  
    Float_t HitposAVRmmU[NMaxHit];         // average of enter and exit position in local ref. frame, U coordinate, in mm
    Float_t HitposAVRmmV[NMaxHit];         // average of enter and exit position in local ref. frame, V coordinate, in mm
    Float_t HitposAVRmmW[NMaxHit];         // average of enter and exit position in local ref. frame, W coordinate, in mm
  
    Float_t HitposAVRmmULadder[NMaxHit];   // always set to -999.0
    Float_t HitposAVRmmVLadder[NMaxHit];   // always set to -999.0
   
    Float_t HitmomMeVX[NMaxHit];           // enter momentum of particle producing the hit, X coordinate, in MeV
    Float_t HitmomMeVY[NMaxHit];           // enter momentum of particle producing the hit, Y coordinate, in MeV
    Float_t HitmomMeVZ[NMaxHit];           // enter momentum of particle producing the hit, Z coordinate, in MeV
   
    Float_t HitthetaLoc[NMaxHit];          // enter theta angle, in deg
    Float_t HitphiLoc[NMaxHit];            // enter phi   angle, in deg
    Float_t HitglobalTime[NMaxHit];        // enter global time, in ns
    Float_t HitGeant4EdepMeV[NMaxHit];     // particles deposited energy by ionization, in MeV
  
    Float_t HitClusterSizeCol[NMaxHit];    // cluster size in column direction => only set if digitization turned on
    Float_t HitClusterSizeRow[NMaxHit];    // cluster size in row    direction => only set if digitization turned on
  
    Float_t HitRecoUmm[NMaxHit];           // reconstructed local hit positio, U coordinate, in mm => only set if digitization turned on
    Float_t HitRecoVmm[NMaxHit];           // reconstructed local hit positio, U coordinate, in mm => only set if digitization turned on
  
    Float_t HitRecoULaddermm[NMaxHit];     // always set to -999.0
    Float_t HitRecoVLaddermm[NMaxHit];     // always set to -999.0
  
    Float_t HitPhiPrincipalAxis[NMaxHit];  // phi angle of cluster principal axis, in deg
    Float_t HitRMSPrincipalAxis[NMaxHit];  // rms of cluster principal axis, in um
    Float_t HitThetaScattering[NMaxHit];   // always set to -1
    
    Int_t   HitNPixels[NMaxHit];           // number of pixels turned on by this hit. If digitization is turned off then it is set to -1
    Int_t   HitFirstPixelIdx[NMaxHit];     // index in pixel block of 1st pixel turned on by this hit
  
    //Definition of Hits Branch 
    Int_t   PixelNb;                       // number of turned on pixels
    Int_t   PixelHitIdx[NMaxPixel];        // index in Hit block of hit turning on this pixel. If pixel from noise, then it is set to -1
    Int_t   PixelGlobalIdx[NMaxPixel];     // pixel global index
    Int_t   PixelColumn[NMaxPixel];        // pixel column
    Int_t   PixelRow[NMaxPixel];           // pixel row
    Float_t PixelAnalogCharge[NMaxPixel];  // pixel collected charge, in electrons (ADC) for digital (analogue) output sensis, i.e. AnalysisMode = 3 (2)
    Int_t   PixelSensorID[NMaxPixel];      // sensor ID
    Float_t PixelUmm[NMaxPixel];           // pixel position, U coordinate, in mm
    Float_t PixelVmm[NMaxPixel];           // pixel position, V coordinate, in mm
    Int_t   PixelStatus[NMaxPixel];        // status variabe. It is always 1 at least the pixel info is lost by saturation, in which case is 0
  
    //Definition of Saturation Branch 
    Int_t   SatuNb;                  // number of lines saturated
    Int_t   SatuLinIdx[NMaxSatu];    // line index
    Int_t   SatuSensorID[NMaxSatu];  // sensor ID
   
    //Definition of EventSize Branch 
    Float_t   EventSizeLadder1; //bytes
    Float_t   EventSizeLadder2; //bytes
    
    
    //Non-sensitive Hit variables
    //Definition of Particles Branch
    Int_t    NonSensitiveParticleNb;                         // number of particles hitting non-sensitive volumes in event
    Int_t    NonSensitiveParticleSensitiveIdx[NMaxParticle]; // index of this particle in Mimosa particle list
    Int_t    NonSensitiveParticlepdgID[NMaxParticle];        // particle pdg ID
    Float_t  NonSensitiveParticleTrkVtX[NMaxParticle];       // particle's origin vertex, X coordinate, in mm
    Float_t  NonSensitiveParticleTrkVtY[NMaxParticle];       // particle's origin vertex, Y coordinate, in mm
    Float_t  NonSensitiveParticleTrkVtZ[NMaxParticle];       // particle's origin vertex, Z coordinate, in mm
    Int_t    NonSensitiveParticleNHits[NMaxParticle];        // particle's number of hits
    Int_t    NonSensitiveParticleFirstHitIdx[NMaxParticle];  // index on hit block of particle's 1st hit
    
    //Definition of Hits Branch
    Int_t   NonSensitiveHitNb;                       // Number of hits in non-sensitive volumes in event
    Int_t   NonSensitiveHitParticleIdx[kMaxHits];    // index in particle block of particle producing this hit
    
    Float_t NonSensitiveHitposINmmX[kMaxHits];       // enter position in global ref. frame, X coordinate, in mm
    Float_t NonSensitiveHitposINmmY[kMaxHits];       // enter position in global ref. frame, Y coordinate, in mm
    Float_t NonSensitiveHitposINmmZ[kMaxHits];       // enter position in global ref. frame, Z coordinate, in mm
    Float_t NonSensitiveHitglobalTimeINns[kMaxHits]; // enter global time, in ns
    Float_t NonSensitiveHitmomINMeVX[kMaxHits];      // enter momentum, X coordinate, in MeV
    Float_t NonSensitiveHitmomINMeVY[kMaxHits];      // enter momentum, Y coordinate, in MeV
    Float_t NonSensitiveHitmomINMeVZ[kMaxHits];      // enter momentum, Z coordinate, in MeV
    Float_t NonSensitiveHitEnergyINMeV[kMaxHits];    // enter Energy, in MeV
    Float_t NonSensitiveHitposOUTmmX[kMaxHits];      // exit  position in global ref. frame, X coordinate, in mm
    Float_t NonSensitiveHitposOUTmmY[kMaxHits];      // exit  position in global ref. frame, Y coordinate, in mm
    Float_t NonSensitiveHitposOUTmmZ[kMaxHits];      // exit  position in global ref. frame, Z coordinate, in mm
    Float_t NonSensitiveHitmomOUTMeVX[kMaxHits];     // exit  momentum, X coordinate, in MeV
    Float_t NonSensitiveHitmomOUTMeVY[kMaxHits];     // exit  momentum, Y coordinate, in MeV
    Float_t NonSensitiveHitmomOUTMeVZ[kMaxHits];     // exit  momentum, Z coordinate, in MeV
    Float_t NonSensitiveHitEnergyOUTMeV[kMaxHits];   // exit  Energy, in MeV
    Float_t NonSensitiveHitGeant4EdepMeV[kMaxHits];  // particles deposited energy by ionization, in MeV
    
    G4double DistFactor;  // distance unit
    G4double AngleFactor; // angle unit
    G4double EPFactor;    // energy-momentum unit
    G4double TimeFactor;  // time unit
    Double_t Edep_correction_factor;
    
    Long64_t  EvtNumber;
    
    bool  NoDigitizationInAllSensitivePlanes;
    
    bool  FillNonSensitiveBranch; // boolean to trigger filling and saving of hits in non-sensitive volumes
    
    //Some function
    void  EqualizeAParticleHits(AParticleHits_t &A, AParticleHits_t B); // Function to equal to AParticleHits_t, i.e. A = B;
    Long64_t  RoundOffInt(double Djentry);                              // Function to roundoff a double to integer
    bool IsAlreadyInPixelList(int Col, int Row, int TheSensorID);       // Check if pixel with in sensor with sensor ID and (col,lin) is already turned on
    double GetATan(double Y, double X);                                 // ATan function
    
  public:
  
    MimosaSimuHistoManager(MimosaSimuSetup* TheSetup);
    ~MimosaSimuHistoManager();
   
    //Some functions
    
    void book();                                                             // book output n-tuple
    void save();                                                             // save output n-tuple 
    void FillMimosaArray(MimosaSimuTrackerHit* aHit, G4int index);           // Fill the hits mimosa (sensitive volume) array
    void FillNonSensitiveArray(MimosaSimuTrajectoryHit* aHit, G4int index);  // Fill the hits in non-sensitive volume array
    
    void FillTreeBlocks(void);                                               // Full up the n-tuple blocks
    void FillUpParticleList(std::vector<Hit_t>  HitsList, std::vector<AParticle_t>&  ListOfParticles); // Fill particle list out of hit list
    void FillUpParticleList_Mimosa(void);       // Fill particle list out of hit list for Mimosa (sensitive volumes)
    void FillUpParticleList_NonSensitive(void); // Fill particle list out of hit list for non-sensitive volumes
    
    
    void   InitChargedParticleList(void);
    bool   IsParticleCharged(int pdgID);
    bool   IsGamma(int pdgID);
    
    int    GetIdxInSensitivePlaneList(int aSensorID);
    
    //Digitization functions
    void InitSensitivePlaneList(void);                                                      // Initialization of sensitive plane list
    void DigitizeEvent(void);                                                               // Digitize all the hits in sensitive planes of the event
    void DigitizeHit(int idx_SD, AParticleHits_t& aHit);                                    // Digitize hit 
    void FillPixelsFromNoise(std::vector<AHitPixels_t> &ListOfPixels);                      // Fills list of pixels from noise
    void FillPixelsFromNoise_Analog( int idx_SD, std::vector<AHitPixels_t> &ListOfPixels);  // Fills pixels from noise of digital  output sensors
    void FillPixelsFromNoise_Digital(int idx_SD, std::vector<AHitPixels_t> &ListOfPixels);  // Fills pixels from noise of analogue output sensors
    
    // n-tuple block filling functions
    void FillMimosaBlock(void);        // Fills up the hits in mimosa (sensitive volumes) block
    void FillNonSensitiveBlock(void);  // Fills up the hits in non-sensitive volumes block
    void FillSaturationBlock(void);    // Fills up the satuation block
    void FillEventSizeBlock(void);     // Fills up event size variables
    void FillNtuple();                 // Fills up the n-tuple

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

