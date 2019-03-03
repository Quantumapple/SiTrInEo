#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


   Long64_t nbytes = 0, nb = 0;

   TH2F *longit = new TH2F("longit","Longitudinal plane; z [mm]; y [mm]",1000,-30,30,1000,-20,20);
   TH2F *transv = new TH2F("transv","Transverse plane ; x [mm]; y [mm]",1000,-30,30,1000,-30,30);

   TGraph *lo1 = new TGraph();
   TGraph *lo2 = new TGraph();
   TGraph *lo3 = new TGraph();
   TGraph *lo4 = new TGraph();

   TGraph *tr1 = new TGraph();
   TGraph *tr2 = new TGraph();
   TGraph *tr3 = new TGraph();
   TGraph *tr4 = new TGraph();

   TFile *output = new TFile("test.root","RECREATE");

   // Set B-field magnitude, Length of B-field
   Float_t Bmag = 0.2; // unit: tesla
   Float_t Bleng = 20*1e-3; // unit: meter
   Float_t pT_true = 2.; // 2 MeV

   Int_t count = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=4; jentry<5;jentry++) {
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

                          Float_t alpha1 = atan((HitposAVRmmV[i+1]-HitposAVRmmV[i])/(HitposINmmZ[i+1]-HitposINmmZ[i]));
                          Float_t alpha2 = atan((HitposAVRmmV[i+3]-HitposAVRmmV[i+2])/(HitposINmmZ[i+3]-HitposINmmZ[i+2]));

                          cout << "             First angle: " << alpha1 << ", second angle: " << alpha2 << endl;
                          cout << "               Angle diff: " << fabs(alpha2-alpha1) << endl;

                          Float_t RecoMom1 = sqrt(pow(HitmomMeVX[i],2)+pow(HitmomMeVY[i],2)+pow(HitmomMeVZ[i],2));
                          Float_t RecoMom2 = sqrt(pow(HitmomMeVX[i+1],2)+pow(HitmomMeVY[i+1],2)+pow(HitmomMeVZ[i+1],2));  
                          Float_t RecoMom3 = sqrt(pow(HitmomMeVX[i+2],2)+pow(HitmomMeVY[i+2],2)+pow(HitmomMeVZ[i+2],2));
                          Float_t RecoMom4 = sqrt(pow(HitmomMeVX[i+3],2)+pow(HitmomMeVY[i+3],2)+pow(HitmomMeVZ[i+3],2));

                          float posU[4] = {HitposAVRmmU[i], HitposAVRmmU[i+1], HitposAVRmmU[i+2], HitposAVRmmU[i+3]};
                          float posX[4] = {HitposINmmX[i], HitposINmmX[i+1], HitposINmmX[i+2], HitposINmmX[i+3]};
                          float posY[4] = {HitposINmmY[i], HitposINmmY[i+1], HitposINmmY[i+2], HitposINmmY[i+3]};

                          lo1->SetPoint(0, HitposINmmZ[i], HitposINmmY[i]); 
                          lo2->SetPoint(0, HitposINmmZ[i+1], HitposINmmY[i+1]); 
                          lo3->SetPoint(0, HitposINmmZ[i+2], HitposINmmY[i+2]); 
                          lo4->SetPoint(0, HitposINmmZ[i+3], HitposINmmY[i+3]); 
                          
                          tr1->SetPoint(0, HitposINmmX[i], HitposINmmY[i]); 
                          tr2->SetPoint(0, HitposINmmX[i+1], HitposINmmY[i+1]); 
                          tr3->SetPoint(0, HitposINmmX[i+2], HitposINmmY[i+2]); 
                          tr4->SetPoint(0, HitposINmmX[i+3], HitposINmmY[i+3]); 

                          cout << "1st-2nd deltaPhi: " << deltaPhi(hit1.Phi(), hit2.Phi()) << ", deltaEta: " << hit1.Eta() - hit2.Eta() << ", deltaR: " << deltaR(hit1.Phi(), hit2.Phi(), hit1.Eta(), hit2.Eta() ) << endl;                   
                          cout << "3rd-4th deltaPhi: " << deltaPhi(hit3.Phi(), hit4.Phi()) << ", deltaEta: " << hit3.Eta() - hit4.Eta() << ", deltaR: " << deltaR(hit3.Phi(), hit4.Phi(), hit3.Eta(), hit4.Eta() ) << endl;                   

                          //cout << "     Hit position x on 1st plane: " << HitposINmmX[i] << ", on 2nd plane: " << HitposINmmX[i+1] << ", on 3rd plane: " << HitposINmmX[i+2] << ", on 4th plane: " << HitposINmmX[i+3] << endl;
                          //cout << "     Hit position y on 1st plane: " << HitposINmmY[i] << ", on 2nd plane: " << HitposINmmY[i+1] << ", on 3rd plane: " << HitposINmmY[i+2] << ", on 4th plane: " << HitposINmmY[i+3] << endl;
                          //cout << "     Hit position U on 1st plane: " << HitposAVRmmU[i] << ", on 2nd plane: " << HitposAVRmmU[i+1] << ", on 3rd plane: " << HitposAVRmmU[i+2] << ", on 4th plane: " << HitposAVRmmU[i+3] << endl;
                          //HXY1->Fill(posX[0], posY[0]);
                          //HXY2->Fill(posX[1], posY[1]);
                          //HXY3->Fill(posX[2], posY[2]);
                          //HXY4->Fill(posX[3], posY[3]);


                      }
                  }
              }
          }

      }
      cout << endl;

      if( flag ) count++;

      /*
      TCanvas *canXY = new TCanvas("canXY","",1000,600);
      HXY1->Draw("");
      HXY1->SetMarkerStyle(20);
      HXY1->SetMarkerColor(9);
      HXY1->GetXaxis()->SetTitle("x (mm)");
      HXY1->GetYaxis()->SetTitle("y (mm)");

      HXY2->Draw("same");
      HXY2->SetMarkerStyle(20);
      HXY2->SetMarkerColor(8);

      HXY3->Draw("same");
      HXY3->SetMarkerStyle(20);
      HXY3->SetMarkerColor(7);

      HXY4->Draw("same");
      HXY4->SetMarkerStyle(20);
      HXY4->SetMarkerColor(6);

      canXY->Write();
      canXY->ResetDrawn();
      */

   } // event loop

   
   lo1->SetMarkerSize(1.5);
   lo1->SetMarkerColor(kBlue);
   lo1->SetMarkerStyle(20);
   lo2->SetMarkerSize(1.5);
   lo2->SetMarkerColor(kRed);
   lo2->SetMarkerStyle(20);
   lo3->SetMarkerSize(1.5);
   lo3->SetMarkerColor(kSpring);
   lo3->SetMarkerStyle(20);
   lo4->SetMarkerSize(1.5);
   lo4->SetMarkerColor(kViolet);
   lo4->SetMarkerStyle(20);
   
   tr1->SetMarkerSize(1.5);
   tr1->SetMarkerColor(kBlue);
   tr1->SetMarkerStyle(20);
   tr2->SetMarkerSize(1.5);
   tr2->SetMarkerColor(kRed);
   tr2->SetMarkerStyle(20);
   tr3->SetMarkerSize(1.5);
   tr3->SetMarkerColor(kSpring);
   tr3->SetMarkerStyle(20);
   tr4->SetMarkerSize(1.5);
   tr4->SetMarkerColor(kViolet);
   tr4->SetMarkerStyle(20);

   TCanvas *c1 = new TCanvas("c1", "", 800,800);
   longit->Draw();
   lo1->Draw("p same");
   lo2->Draw("p same");
   lo3->Draw("p same");
   lo4->Draw("p same");

   TCanvas *c2 = new TCanvas("c2", "", 800, 800);
   transv->Draw();
   tr1->Draw("p same");
   tr2->Draw("p same");
   tr3->Draw("p same");
   tr4->Draw("p same");

   c1->Print("LongitGood.png");
   c2->Print("TransvGood.png");

   delete c1;
   delete c2;
   
   //cout << endl;
   //cout << "Total event : "<< nentries << endl;
   //cout << "   Number of events at least one track measured with 4 hits on each pixel layer: " << count << endl;
   //cout << "     Signal efficiency: " << (float)count/3000. << endl;

   output->Close();

}
