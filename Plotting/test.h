//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 26 21:54:32 2019 by ROOT version 6.12/06
// from TTree ntp1/Hits
// found on file: runXXX.cfg.BT.sitrineo_d5_l5_s10_B02_L20.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ParticleNb;
   Int_t           ParticlepdgID[6];   //[ParticleNb]
   Int_t           ParticleBKGType[6];   //[ParticleNb]
   Float_t         ParticleTrkVtX[6];   //[ParticleNb]
   Float_t         ParticleTrkVtY[6];   //[ParticleNb]
   Float_t         ParticleTrkVtZ[6];   //[ParticleNb]
   Int_t           ParticleNHits[6];   //[ParticleNb]
   Int_t           ParticleFirstHitIdx[6];   //[ParticleNb]
   Int_t           HitNb;
   Int_t           HitParticleIdx[23];   //[HitNb]
   Int_t           HitsensorID[23];   //[HitNb]
   Int_t           HitladderID[23];   //[HitNb]
   Int_t           HitmoduleID[23];   //[HitNb]
   Float_t         HitposINmmX[23];   //[HitNb]
   Float_t         HitposINmmY[23];   //[HitNb]
   Float_t         HitposINmmZ[23];   //[HitNb]
   Float_t         HitposAVRmmU[23];   //[HitNb]
   Float_t         HitposAVRmmV[23];   //[HitNb]
   Float_t         HitposAVRmmW[23];   //[HitNb]
   Float_t         HitposAVRmmULadder[23];   //[HitNb]
   Float_t         HitposAVRmmVLadder[23];   //[HitNb]
   Float_t         HitmomMeVX[23];   //[HitNb]
   Float_t         HitmomMeVY[23];   //[HitNb]
   Float_t         HitmomMeVZ[23];   //[HitNb]
   Float_t         HitthetaLoc[23];   //[HitNb]
   Float_t         HitphiLoc[23];   //[HitNb]
   Float_t         HitglobalTime[23];   //[HitNb]
   Float_t         HitGeant4EdepMeV[23];   //[HitNb]
   Float_t         HitClusterSizeCol[23];   //[HitNb]
   Float_t         HitClusterSizeRow[23];   //[HitNb]
   Float_t         HitRecoUmm[23];   //[HitNb]
   Float_t         HitRecoVmm[23];   //[HitNb]
   Float_t         HitRecoULaddermm[23];   //[HitNb]
   Float_t         HitRecoVLaddermm[23];   //[HitNb]
   Float_t         HitPhiPrincipalAxis[23];   //[HitNb]
   Float_t         HitRMSPrincipalAxis[23];   //[HitNb]
   Float_t         HitThetaScattering[23];   //[HitNb]
   Int_t           HitNPixels[23];   //[HitNb]
   Int_t           HitFirstPixelIdx[23];   //[HitNb]
   Int_t           PixelNb;
   Int_t           PixelHitIdx[336];   //[PixelNb]
   Int_t           PixelGlobalIdx[336];   //[PixelNb]
   Int_t           PixelColumn[336];   //[PixelNb]
   Int_t           PixelRow[336];   //[PixelNb]
   Float_t         PixelAnalogCharge[336];   //[PixelNb]
   Int_t           PixelSensorID[336];   //[PixelNb]
   Float_t         PixelUmm[336];   //[PixelNb]
   Float_t         PixelVmm[336];   //[PixelNb]
   Int_t           PixelStatus[336];   //[PixelNb]
   Int_t           SaturationInfoNb;
   Int_t           SaturationInfoLinIdx[1];   //[SaturationInfoNb]
   Int_t           SaturationInfoSensorID[1];   //[SaturationInfoNb]
   Float_t         EventSizeLadder1;
   Float_t         EventSizeLadder2;
   Int_t           NonSensitiveParticleNb;
   Int_t           NonSensitiveParticleSensitiveIdx[13];   //[NonSensitiveParticleNb]
   Int_t           NonSensitiveParticlepdgID[13];   //[NonSensitiveParticleNb]
   Float_t         NonSensitiveParticleTrkVtX[13];   //[NonSensitiveParticleNb]
   Float_t         NonSensitiveParticleTrkVtY[13];   //[NonSensitiveParticleNb]
   Float_t         NonSensitiveParticleTrkVtZ[13];   //[NonSensitiveParticleNb]
   Int_t           NonSensitiveParticleNHits[13];   //[NonSensitiveParticleNb]
   Int_t           NonSensitiveParticleFirstHitIdx[13];   //[NonSensitiveParticleNb]
   Int_t           NonSensitiveHitNb;
   Int_t           NonSensitiveHitParticleIdx[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitposINmmX[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitposINmmY[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitposINmmZ[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitglobalTimeINns[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitmomINMeVX[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitmomINMeVY[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitmomINMeVZ[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitEnergyINMeV[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitposOUTmmX[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitposOUTmmY[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitposOUTmmZ[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitmomOUTMeVX[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitmomOUTMeVY[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitmomOUTMeVZ[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitEnergyOUTMeV[668];   //[NonSensitiveHitNb]
   Float_t         NonSensitiveHitGeant4EdepMeV[668];   //[NonSensitiveHitNb]

   // List of branches
   TBranch        *b_ParticleNb;   //!
   TBranch        *b_ParticlepdgID;   //!
   TBranch        *b_ParticleBKGType;   //!
   TBranch        *b_ParticleTrkVtX;   //!
   TBranch        *b_ParticleTrkVtY;   //!
   TBranch        *b_ParticleTrkVtZ;   //!
   TBranch        *b_ParticleNHits;   //!
   TBranch        *b_ParticleFirstHitIdx;   //!
   TBranch        *b_HitNb;   //!
   TBranch        *b_HitParticleIdx;   //!
   TBranch        *b_HitsensorID;   //!
   TBranch        *b_HitladderID;   //!
   TBranch        *b_HitmoduleID;   //!
   TBranch        *b_HitposINmmX;   //!
   TBranch        *b_HitposINmmY;   //!
   TBranch        *b_HitposINmmZ;   //!
   TBranch        *b_HitposAVRmmU;   //!
   TBranch        *b_HitposAVRmmV;   //!
   TBranch        *b_HitposAVRmmW;   //!
   TBranch        *b_HitposAVRmmULadder;   //!
   TBranch        *b_HitposAVRmmVLadder;   //!
   TBranch        *b_HitmomMeVX;   //!
   TBranch        *b_HitmomMeVY;   //!
   TBranch        *b_HitmomMeVZ;   //!
   TBranch        *b_HitthetaLoc;   //!
   TBranch        *b_HitphiLoc;   //!
   TBranch        *b_HitglobalTime;   //!
   TBranch        *b_HitGeant4EdepMeV;   //!
   TBranch        *b_HitClusterSizeCol;   //!
   TBranch        *b_HitClusterSizeRow;   //!
   TBranch        *b_HitRecoUmm;   //!
   TBranch        *b_HitRecoVmm;   //!
   TBranch        *b_HitRecoULaddermm;   //!
   TBranch        *b_HitRecoVLaddermm;   //!
   TBranch        *b_HitPhiPrincipalAxis;   //!
   TBranch        *b_HitRMSPrincipalAxis;   //!
   TBranch        *b_HitThetaScattering;   //!
   TBranch        *b_HitNPixels;   //!
   TBranch        *b_HitFirstPixelIdx;   //!
   TBranch        *b_PixelNb;   //!
   TBranch        *b_PixelHitIdx;   //!
   TBranch        *b_PixelGlobalIdx;   //!
   TBranch        *b_PixelColumn;   //!
   TBranch        *b_PixelRow;   //!
   TBranch        *b_PixelAnalogCharge;   //!
   TBranch        *b_PixelSensorID;   //!
   TBranch        *b_PixelUmm;   //!
   TBranch        *b_PixelVmm;   //!
   TBranch        *b_PixelStatus;   //!
   TBranch        *b_SaturationInfoNb;   //!
   TBranch        *b_SaturationInfoLinIdx;   //!
   TBranch        *b_SaturationInfoSensorID;   //!
   TBranch        *b_EventSizeLadder1;   //!
   TBranch        *b_EventSizeLadder2;   //!
   TBranch        *b_NonSensitiveParticleNb;   //!
   TBranch        *b_NonSensitiveParticleSensitiveIdx;   //!
   TBranch        *b_NonSensitiveParticlepdgID;   //!
   TBranch        *b_NonSensitiveParticleTrkVtX;   //!
   TBranch        *b_NonSensitiveParticleTrkVtY;   //!
   TBranch        *b_NonSensitiveParticleTrkVtZ;   //!
   TBranch        *b_NonSensitiveParticleNHits;   //!
   TBranch        *b_NonSensitiveParticleFirstHitIdx;   //!
   TBranch        *b_NonSensitiveHitNb;   //!
   TBranch        *b_NonSensitiveHitParticleIdx;   //!
   TBranch        *b_NonSensitiveHitposINmmX;   //!
   TBranch        *b_NonSensitiveHitposINmmY;   //!
   TBranch        *b_NonSensitiveHitposINmmZ;   //!
   TBranch        *b_NonSensitiveHitglobalTimeINns;   //!
   TBranch        *b_NonSensitiveHitmomINMeVX;   //!
   TBranch        *b_NonSensitiveHitmomINMeVY;   //!
   TBranch        *b_NonSensitiveHitmomINMeVZ;   //!
   TBranch        *b_NonSensitiveHitEnergyINMeV;   //!
   TBranch        *b_NonSensitiveHitposOUTmmX;   //!
   TBranch        *b_NonSensitiveHitposOUTmmY;   //!
   TBranch        *b_NonSensitiveHitposOUTmmZ;   //!
   TBranch        *b_NonSensitiveHitmomOUTMeVX;   //!
   TBranch        *b_NonSensitiveHitmomOUTMeVY;   //!
   TBranch        *b_NonSensitiveHitmomOUTMeVZ;   //!
   TBranch        *b_NonSensitiveHitEnergyOUTMeV;   //!
   TBranch        *b_NonSensitiveHitGeant4EdepMeV;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   inline float deltaPhi(float phi1, float phi2)
   {
       Float_t tmpDp = phi1 - phi2;
       if( tmpDp >= TMath::Pi() ) tmpDp -= 2.*TMath::Pi();
       if( tmpDp < -TMath::Pi() ) tmpDp += 2.*TMath::Pi();
       return tmpDp;
   }


   inline float deltaR(float phi1, float phi2, float eta1, float eta2)
   {
       return sqrt(pow(eta2-eta1,2)+pow(deltaPhi(phi1, phi2),2));
   }

   inline float SqRt(float var1, float var2)
   {
       return sqrt(pow(var1,2)+pow(var2,2));
   }

};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("runXXX.cfg.BT.sitrineo_d5_l5_s10_B02_L20.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("runXXX.cfg.BT.sitrineo_d5_l5_s10_B02_L20.root");
      }
      f->GetObject("ntp1",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ParticleNb", &ParticleNb, &b_ParticleNb);
   fChain->SetBranchAddress("ParticlepdgID", ParticlepdgID, &b_ParticlepdgID);
   fChain->SetBranchAddress("ParticleBKGType", ParticleBKGType, &b_ParticleBKGType);
   fChain->SetBranchAddress("ParticleTrkVtX", ParticleTrkVtX, &b_ParticleTrkVtX);
   fChain->SetBranchAddress("ParticleTrkVtY", ParticleTrkVtY, &b_ParticleTrkVtY);
   fChain->SetBranchAddress("ParticleTrkVtZ", ParticleTrkVtZ, &b_ParticleTrkVtZ);
   fChain->SetBranchAddress("ParticleNHits", ParticleNHits, &b_ParticleNHits);
   fChain->SetBranchAddress("ParticleFirstHitIdx", ParticleFirstHitIdx, &b_ParticleFirstHitIdx);
   fChain->SetBranchAddress("HitNb", &HitNb, &b_HitNb);
   fChain->SetBranchAddress("HitParticleIdx", HitParticleIdx, &b_HitParticleIdx);
   fChain->SetBranchAddress("HitsensorID", HitsensorID, &b_HitsensorID);
   fChain->SetBranchAddress("HitladderID", HitladderID, &b_HitladderID);
   fChain->SetBranchAddress("HitmoduleID", HitmoduleID, &b_HitmoduleID);
   fChain->SetBranchAddress("HitposINmmX", HitposINmmX, &b_HitposINmmX);
   fChain->SetBranchAddress("HitposINmmY", HitposINmmY, &b_HitposINmmY);
   fChain->SetBranchAddress("HitposINmmZ", HitposINmmZ, &b_HitposINmmZ);
   fChain->SetBranchAddress("HitposAVRmmU", HitposAVRmmU, &b_HitposAVRmmU);
   fChain->SetBranchAddress("HitposAVRmmV", HitposAVRmmV, &b_HitposAVRmmV);
   fChain->SetBranchAddress("HitposAVRmmW", HitposAVRmmW, &b_HitposAVRmmW);
   fChain->SetBranchAddress("HitposAVRmmULadder", HitposAVRmmULadder, &b_HitposAVRmmULadder);
   fChain->SetBranchAddress("HitposAVRmmVLadder", HitposAVRmmVLadder, &b_HitposAVRmmVLadder);
   fChain->SetBranchAddress("HitmomMeVX", HitmomMeVX, &b_HitmomMeVX);
   fChain->SetBranchAddress("HitmomMeVY", HitmomMeVY, &b_HitmomMeVY);
   fChain->SetBranchAddress("HitmomMeVZ", HitmomMeVZ, &b_HitmomMeVZ);
   fChain->SetBranchAddress("HitthetaLoc", HitthetaLoc, &b_HitthetaLoc);
   fChain->SetBranchAddress("HitphiLoc", HitphiLoc, &b_HitphiLoc);
   fChain->SetBranchAddress("HitglobalTime", HitglobalTime, &b_HitglobalTime);
   fChain->SetBranchAddress("HitGeant4EdepMeV", HitGeant4EdepMeV, &b_HitGeant4EdepMeV);
   fChain->SetBranchAddress("HitClusterSizeCol", HitClusterSizeCol, &b_HitClusterSizeCol);
   fChain->SetBranchAddress("HitClusterSizeRow", HitClusterSizeRow, &b_HitClusterSizeRow);
   fChain->SetBranchAddress("HitRecoUmm", HitRecoUmm, &b_HitRecoUmm);
   fChain->SetBranchAddress("HitRecoVmm", HitRecoVmm, &b_HitRecoVmm);
   fChain->SetBranchAddress("HitRecoULaddermm", HitRecoULaddermm, &b_HitRecoULaddermm);
   fChain->SetBranchAddress("HitRecoVLaddermm", HitRecoVLaddermm, &b_HitRecoVLaddermm);
   fChain->SetBranchAddress("HitPhiPrincipalAxis", HitPhiPrincipalAxis, &b_HitPhiPrincipalAxis);
   fChain->SetBranchAddress("HitRMSPrincipalAxis", HitRMSPrincipalAxis, &b_HitRMSPrincipalAxis);
   fChain->SetBranchAddress("HitThetaScattering", HitThetaScattering, &b_HitThetaScattering);
   fChain->SetBranchAddress("HitNPixels", HitNPixels, &b_HitNPixels);
   fChain->SetBranchAddress("HitFirstPixelIdx", HitFirstPixelIdx, &b_HitFirstPixelIdx);
   fChain->SetBranchAddress("PixelNb", &PixelNb, &b_PixelNb);
   fChain->SetBranchAddress("PixelHitIdx", PixelHitIdx, &b_PixelHitIdx);
   fChain->SetBranchAddress("PixelGlobalIdx", PixelGlobalIdx, &b_PixelGlobalIdx);
   fChain->SetBranchAddress("PixelColumn", PixelColumn, &b_PixelColumn);
   fChain->SetBranchAddress("PixelRow", PixelRow, &b_PixelRow);
   fChain->SetBranchAddress("PixelAnalogCharge", PixelAnalogCharge, &b_PixelAnalogCharge);
   fChain->SetBranchAddress("PixelSensorID", PixelSensorID, &b_PixelSensorID);
   fChain->SetBranchAddress("PixelUmm", PixelUmm, &b_PixelUmm);
   fChain->SetBranchAddress("PixelVmm", PixelVmm, &b_PixelVmm);
   fChain->SetBranchAddress("PixelStatus", PixelStatus, &b_PixelStatus);
   fChain->SetBranchAddress("SaturationInfoNb", &SaturationInfoNb, &b_SaturationInfoNb);
   fChain->SetBranchAddress("SaturationInfoLinIdx", &SaturationInfoLinIdx, &b_SaturationInfoLinIdx);
   fChain->SetBranchAddress("SaturationInfoSensorID", &SaturationInfoSensorID, &b_SaturationInfoSensorID);
   fChain->SetBranchAddress("EventSizeLadder1", &EventSizeLadder1, &b_EventSizeLadder1);
   fChain->SetBranchAddress("EventSizeLadder2", &EventSizeLadder2, &b_EventSizeLadder2);
   fChain->SetBranchAddress("NonSensitiveParticleNb", &NonSensitiveParticleNb, &b_NonSensitiveParticleNb);
   fChain->SetBranchAddress("NonSensitiveParticleSensitiveIdx", NonSensitiveParticleSensitiveIdx, &b_NonSensitiveParticleSensitiveIdx);
   fChain->SetBranchAddress("NonSensitiveParticlepdgID", NonSensitiveParticlepdgID, &b_NonSensitiveParticlepdgID);
   fChain->SetBranchAddress("NonSensitiveParticleTrkVtX", NonSensitiveParticleTrkVtX, &b_NonSensitiveParticleTrkVtX);
   fChain->SetBranchAddress("NonSensitiveParticleTrkVtY", NonSensitiveParticleTrkVtY, &b_NonSensitiveParticleTrkVtY);
   fChain->SetBranchAddress("NonSensitiveParticleTrkVtZ", NonSensitiveParticleTrkVtZ, &b_NonSensitiveParticleTrkVtZ);
   fChain->SetBranchAddress("NonSensitiveParticleNHits", NonSensitiveParticleNHits, &b_NonSensitiveParticleNHits);
   fChain->SetBranchAddress("NonSensitiveParticleFirstHitIdx", NonSensitiveParticleFirstHitIdx, &b_NonSensitiveParticleFirstHitIdx);
   fChain->SetBranchAddress("NonSensitiveHitNb", &NonSensitiveHitNb, &b_NonSensitiveHitNb);
   fChain->SetBranchAddress("NonSensitiveHitParticleIdx", NonSensitiveHitParticleIdx, &b_NonSensitiveHitParticleIdx);
   fChain->SetBranchAddress("NonSensitiveHitposINmmX", NonSensitiveHitposINmmX, &b_NonSensitiveHitposINmmX);
   fChain->SetBranchAddress("NonSensitiveHitposINmmY", NonSensitiveHitposINmmY, &b_NonSensitiveHitposINmmY);
   fChain->SetBranchAddress("NonSensitiveHitposINmmZ", NonSensitiveHitposINmmZ, &b_NonSensitiveHitposINmmZ);
   fChain->SetBranchAddress("NonSensitiveHitglobalTimeINns", NonSensitiveHitglobalTimeINns, &b_NonSensitiveHitglobalTimeINns);
   fChain->SetBranchAddress("NonSensitiveHitmomINMeVX", NonSensitiveHitmomINMeVX, &b_NonSensitiveHitmomINMeVX);
   fChain->SetBranchAddress("NonSensitiveHitmomINMeVY", NonSensitiveHitmomINMeVY, &b_NonSensitiveHitmomINMeVY);
   fChain->SetBranchAddress("NonSensitiveHitmomINMeVZ", NonSensitiveHitmomINMeVZ, &b_NonSensitiveHitmomINMeVZ);
   fChain->SetBranchAddress("NonSensitiveHitEnergyINMeV", NonSensitiveHitEnergyINMeV, &b_NonSensitiveHitEnergyINMeV);
   fChain->SetBranchAddress("NonSensitiveHitposOUTmmX", NonSensitiveHitposOUTmmX, &b_NonSensitiveHitposOUTmmX);
   fChain->SetBranchAddress("NonSensitiveHitposOUTmmY", NonSensitiveHitposOUTmmY, &b_NonSensitiveHitposOUTmmY);
   fChain->SetBranchAddress("NonSensitiveHitposOUTmmZ", NonSensitiveHitposOUTmmZ, &b_NonSensitiveHitposOUTmmZ);
   fChain->SetBranchAddress("NonSensitiveHitmomOUTMeVX", NonSensitiveHitmomOUTMeVX, &b_NonSensitiveHitmomOUTMeVX);
   fChain->SetBranchAddress("NonSensitiveHitmomOUTMeVY", NonSensitiveHitmomOUTMeVY, &b_NonSensitiveHitmomOUTMeVY);
   fChain->SetBranchAddress("NonSensitiveHitmomOUTMeVZ", NonSensitiveHitmomOUTMeVZ, &b_NonSensitiveHitmomOUTMeVZ);
   fChain->SetBranchAddress("NonSensitiveHitEnergyOUTMeV", NonSensitiveHitEnergyOUTMeV, &b_NonSensitiveHitEnergyOUTMeV);
   fChain->SetBranchAddress("NonSensitiveHitGeant4EdepMeV", NonSensitiveHitGeant4EdepMeV, &b_NonSensitiveHitGeant4EdepMeV);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
