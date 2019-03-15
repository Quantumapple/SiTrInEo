#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

// Check MimosaSimuHistoManager.hh to understand the branch

TH1F *hist_Style(TH1F *h, TString xtitle, TString ytitle)
{
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    TH1F *hnew = (TH1F*)h->Clone("new");
    return hnew;
}


void test::Loop()
{
   if (fChain == 0) return;

   gStyle->SetOptStat(1110);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TH1F *dR32 = new TH1F("dR32", "#Delta R distribution between 1st and 2nd pixel hits", 200, 0., 5.);
   TH1F *dR21 = new TH1F("dR21", "#Delta R distribution between 2nd and 3rd pixel hits", 200, 0., 10.);
   TH1F *dR10 = new TH1F("dR10", "#Delta R distribution between 3rd and 4th pixel hits", 200, 0., 5.);

   TH1F *diffU = new TH1F("diffU", "Difference AVRUmm and RecoUmm", 500, -0.05, 0.05);
   TH1F *diffV = new TH1F("diffV", "Difference AVRVmm and RecoVmm", 500, -0.05, 0.05);

   TH1F *MomPix1 = new TH1F("MomPix1", "Momentum distribution of 1st pixel hits", 100, 0, 2.5);
   TH1F *MomPix2 = new TH1F("MomPix2", "Momentum distribution of 2nd pixel hits", 100, 0, 2.5);
   TH1F *MomPix3 = new TH1F("MomPix3", "Momentum distribution of 3rd pixel hits", 100, 0, 2.5);
   TH1F *MomPix4 = new TH1F("MomPix4", "Momentum distribution of 4th pixel hits", 100, 0, 2.5);

   TH1F *reCoPt = new TH1F("reCoPt", "Reconstructed momentum from the track", 100, 0, 10.);
   //TH1F *resiPt = new TH1F("resiPt", "Residual of pT", 100, 0, 10.);
   TH2F *resiPt = new TH2F("resiPt", "Residual of pT", 100, 0., 2.5, 1000, 0., 10.);

   TH1F *ang1 = new TH1F("ang1", "Angle from 1st and 2nd pixels", 200, -2, 2);
   TH1F *ang2 = new TH1F("ang2", "Angle from 3rd and 4th pixels", 200, -2, 2);
   TH1F *diffAng = new TH1F("diffAng", "Angle difference of two angles", 200, -2, 2);

   // Set B-field magnitude, Length of B-field
   Float_t Bmag = 0.2; // unit: tesla
   Float_t Bleng = 20*1e-3; // unit: meter
   Float_t pT_true = 2.; // 2 MeV

   Int_t count = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<50;jentry++) {
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
          diffU->Fill(HitposAVRmmU[i]-HitRecoUmm[i]);
          diffV->Fill(HitposAVRmmV[i]-HitRecoVmm[i]);

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

                          Float_t alpha1 = atan((HitposAVRmmV[i+1]-HitposAVRmmV[i])/(HitposINmmZ[i+1]-HitposINmmZ[i]));
                          Float_t alpha2 = atan((HitposAVRmmV[i+3]-HitposAVRmmV[i+2])/(HitposINmmZ[i+3]-HitposINmmZ[i+2]));

                          cout << "             First angle: " << alpha1 << ", second angle: " << alpha2 << endl;
                          cout << "               Angle diff: " << fabs(alpha2-alpha1) << endl;
                          ang1->Fill(alpha1); ang2->Fill(alpha2);
                          diffAng->Fill(alpha2-alpha1);

                          Float_t RecoMom1 = sqrt(pow(HitmomMeVX[i],2)+pow(HitmomMeVY[i],2)+pow(HitmomMeVZ[i],2));
                          Float_t RecoMom2 = sqrt(pow(HitmomMeVX[i+1],2)+pow(HitmomMeVY[i+1],2)+pow(HitmomMeVZ[i+1],2));  
                          Float_t RecoMom3 = sqrt(pow(HitmomMeVX[i+2],2)+pow(HitmomMeVY[i+2],2)+pow(HitmomMeVZ[i+2],2));
                          Float_t RecoMom4 = sqrt(pow(HitmomMeVX[i+3],2)+pow(HitmomMeVY[i+3],2)+pow(HitmomMeVZ[i+3],2));

                          MomPix1->Fill(RecoMom1); 
                          MomPix2->Fill(RecoMom2); 
                          MomPix3->Fill(RecoMom3); 
                          MomPix4->Fill(RecoMom4);

                          Float_t RecoPt = (0.3*Bmag*Bleng)/fabs(alpha2-alpha1)*1e+3;
                          reCoPt->Fill(RecoPt);

                          Float_t diffPt = fabs(RecoPt - pT_true);
                          resiPt->Fill(pT_true, diffPt/pT_true);

                      }
                  }
              }
          }

      }
      cout << endl;

      if( flag ) count++;
   
   } // event loop

   diffU->GetXaxis()->SetTitle("#Delta U [mm]");
   diffV->GetXaxis()->SetTitle("#Delta V [mm]");

   //cout << endl;
   //cout << "Total event: 3000" << endl;
   //cout << "   Number of events at least one track measured with 4 hits on each pixel layer: " << count << endl;
   //cout << "     Signal efficiency: " << (float)count/3000. << endl;

   Int_t total = dR32->Integral();
   cout << "Check the signal efficiency based on Delta R" << endl;
   cout << "  From 1st - 2nd Delta R histo" << endl;
   cout << "    When Delta R = 0.8, left: " << (float)dR32->Integral(1,16)/total << ", right: " << 1. - (float)dR32->Integral(1,16)/total << endl;
   cout << "    When Delta R = 0.9, left: " << (float)dR32->Integral(1,18)/total << ", right: " << 1. - (float)dR32->Integral(1,18)/total << endl;
   cout << "    When Delta R = 1.0, left: " << (float)dR32->Integral(1,20)/total << ", right: " << 1. - (float)dR32->Integral(1,20)/total << endl;
   cout << "    When Delta R = 1.1, left: " << (float)dR32->Integral(1,22)/total << ", right: " << 1. - (float)dR32->Integral(1,22)/total << endl;
   cout << "    When Delta R = 1.2, left: " << (float)dR32->Integral(1,24)/total << ", right: " << 1. - (float)dR32->Integral(1,24)/total << endl;
   cout << "    When Delta R = 1.3, left: " << (float)dR32->Integral(1,26)/total << ", right: " << 1. - (float)dR32->Integral(1,26)/total << endl;

   cout << "  From 3rd - 4th Delta R histo" << endl;
   cout << "    When Delta R = 0.5, left: " << (float)dR10->Integral(1,10)/total << ", right: " << 1. - (float)dR10->Integral(1,10)/total << endl;
   cout << "    When Delta R = 0.6, left: " << (float)dR10->Integral(1,12)/total << ", right: " << 1. - (float)dR10->Integral(1,12)/total << endl;
   cout << "    When Delta R = 0.7, left: " << (float)dR10->Integral(1,14)/total << ", right: " << 1. - (float)dR10->Integral(1,14)/total << endl;
   cout << "    When Delta R = 0.8, left: " << (float)dR10->Integral(1,16)/total << ", right: " << 1. - (float)dR10->Integral(1,16)/total << endl;
   cout << "    When Delta R = 0.9, left: " << (float)dR10->Integral(1,18)/total << ", right: " << 1. - (float)dR10->Integral(1,18)/total << endl;
   cout << "    When Delta R = 1.0, left: " << (float)dR10->Integral(1,20)/total << ", right: " << 1. - (float)dR10->Integral(1,20)/total << endl;
   cout << "    When Delta R = 1.1, left: " << (float)dR10->Integral(1,22)/total << ", right: " << 1. - (float)dR10->Integral(1,22)/total << endl;
   cout << "    When Delta R = 1.2, left: " << (float)dR10->Integral(1,24)/total << ", right: " << 1. - (float)dR10->Integral(1,24)/total << endl;
   cout << "    When Delta R = 1.3, left: " << (float)dR10->Integral(1,26)/total << ", right: " << 1. - (float)dR10->Integral(1,26)/total << endl;

   TString Ytitle = "Events";
   hist_Style(dR32, "#Delta R", Ytitle);
   hist_Style(dR10, "#Delta R", Ytitle);
  
   hist_Style(MomPix1, "p_{T}", Ytitle);
   hist_Style(MomPix2, "p_{T}", Ytitle);
   hist_Style(MomPix3, "p_{T}", Ytitle);
   hist_Style(MomPix4, "p_{T}", Ytitle);
   hist_Style(reCoPt, "p^{reco}_{T}", Ytitle);

   hist_Style(ang1, "#alpha_{1}(1st, 2nd)", Ytitle);
   hist_Style(ang2, "#alpha_{2}(3rd, 4th)", Ytitle);
   hist_Style(diffAng, "#Delta#alpha = #alpha_{2} - #alpha_{1}", Ytitle);


   TCanvas *c1 = new TCanvas("c1","",600,600);
   dR32->Draw();
   c1->Print("dR32.png");
   dR10->Draw();
   c1->Print("dR10.png");

   MomPix1->SetLineColor(kBlue);
   MomPix2->SetLineColor(kRed);
   MomPix3->SetLineColor(kSpring);
   MomPix4->SetLineColor(kViolet);

   MomPix1->Draw();
   //c1->Print("InitPt1.png");
   MomPix2->Draw("same");
   //c1->Print("InitPt2.png");
   MomPix3->Draw("same");
   //c1->Print("InitPt3.png");
   MomPix4->Draw("same");

   TLegend *lgd1 = new TLegend(0.2, 0.4, 0.65, 0.79);
   lgd1->AddEntry(MomPix1, "Initial p_{T} from 1st pixel", "l");
   lgd1->AddEntry(MomPix2, "Initial p_{T} from 2nd pixel", "l");
   lgd1->AddEntry(MomPix3, "Initial p_{T} from 3rd pixel", "l");
   lgd1->AddEntry(MomPix4, "Initial p_{T} from 4th pixel", "l");
   lgd1->Draw();
   c1->Print("InitPt.png");

   reCoPt->Draw();
   c1->Print("RecoPt.png");

   ang1->Draw();
   c1->Print("ang1.png");
   ang2->Draw();
   c1->Print("ang2.png");
   diffAng->Draw();
   c1->Print("diffAng.png");

   TFile *output = new TFile("output.root","RECREATE");
   dR32->Write();
   dR21->Write();
   dR10->Write();

   diffU->Write();
   diffV->Write();
   
   MomPix1->Write(); 
   MomPix2->Write(); 
   MomPix3->Write(); 
   MomPix4->Write(); 

   reCoPt->Write();

   resiPt->Write();

   ang1->Write();
   ang2->Write();
   diffAng->Write();
}
