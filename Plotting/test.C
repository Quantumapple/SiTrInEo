#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

// Check MimosaSimuHistoManager.hh to understand the branch

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH1F *dR32 = new TH1F("dR32", "#Delta R distribution between 1st and 2nd pixel hits", 200, 0., 10.);
   TH1F *dR21 = new TH1F("dR21", "#Delta R distribution between 2nd and 3rd pixel hits", 200, 0., 10.);
   TH1F *dR10 = new TH1F("dR10", "#Delta R distribution between 3rd and 4th pixel hits", 200, 0., 10.);

   TH1F *ang1 = new TH1F("ang1", "Angle from 1st and 2nd pixel hits", 100, 0., 3.15);
   TH1F *ang2 = new TH1F("ang2", "Angle from 3rd and 4th pixel hits", 100, 0., 3.15);

   Int_t count = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      cout << jentry << " event" << endl;
      //cout << "    How many particles? " << ParticleNb << endl; 
      //cout << endl;
      //for(Int_t i = 0; i < ParticleNb; i++)
      //{
      //}

      Bool_t flag = false;
      for(Int_t i = 0; i < HitNb; i++)
      {
          //cout << "     Which sensor? " << HitsensorID[i] << endl;
          if( HitsensorID[i] == 3)
          {
              cout << "  Which sensor? " << HitsensorID[i] << endl;
              TVector3 hit1;
              hit1.SetXYZ(HitposINmmX[i], HitposINmmY[i], HitposINmmZ[i]);
              Float_t phi1 = hit1.Phi();
              Float_t eta1 = hit1.Eta();
              if( HitsensorID[i+1] == 2 )
              {
                  cout << "    Which sensor? " << HitsensorID[i+1] << endl;
                  TVector3 hit2;
                  hit2.SetXYZ(HitposINmmX[i+1], HitposINmmY[i+1], HitposINmmZ[i+1]);
                  Float_t phi2 = hit2.Phi();
                  Float_t eta2 = hit2.Eta();
                  if( HitsensorID[i+2] == 1 )
                  {
                      cout << "      Which sensor? " << HitsensorID[i+2] << endl;
                      TVector3 hit3;
                      hit3.SetXYZ(HitposINmmX[i+2], HitposINmmY[i+2], HitposINmmZ[i+2]);
                      Float_t phi3 = hit3.Phi();
                      Float_t eta3 = hit3.Eta();
                      if( HitsensorID[i+3] == 0 )
                      {
                          TVector3 hit4;
                          hit4.SetXYZ(HitposINmmX[i+3], HitposINmmY[i+3], HitposINmmZ[i+3]);
                          Float_t phi4 = hit4.Phi();
                          Float_t eta4 = hit4.Eta();
                          cout << "        Which sensor? " << HitsensorID[i+3] << endl;
                          cout << "           This event has a track which has 4 hits on each pixel layer" << endl;
                          //count++;
                          flag = true;

                          dR32->Fill(deltaR(phi1, phi2, eta1, eta2));
                          dR21->Fill(deltaR(phi2, phi3, eta2, eta3));
                          dR10->Fill(deltaR(phi3, phi4, eta3, eta4));

                          Float_t alpha1 = 

                      }
                  }
              }
          }

      }
      cout << endl;

      if( flag ) count++;
   
   } // event loop

   cout << endl;
   cout << "Total event: 3000" << endl;
   cout << "   Number of events at least one track measured with 4 hits on each pixel layer: " << count << endl;
   cout << "     Signal efficiency: " << (float)count/3000. << endl;

   TFile *output = new TFile("output.root","RECREATE");
   dR32->Write();
   dR21->Write();
   dR10->Write();
}
