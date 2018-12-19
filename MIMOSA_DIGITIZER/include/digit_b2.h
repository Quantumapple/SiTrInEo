#ifndef DIGIT_B2_H
#define DIGIT_B2_H
///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                           //
//     /home/abesson/SAVE/ILCSOFT/DIGITISEUR_BELLEII                                         //
//                                                                                           //
//       {
//        gROOT->ProcessLine(".L digit_b2.cxx+");
//        DIGIT_B2 myDIGIT_B2("myname","mytitle");
//       }
//                                                                                           //
//                                                                                           //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TList.h>
#include <TGraph.h>
#include <TCanvas.h>
#include "Riostream.h"
#include "vector"
#include <TRandom3.h>

// ROOT classes
#include "TStopwatch.h"
#include "TString.h"
#include "TObject.h"
#include "TVector.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"

using namespace std;

//Int_t GlobalSeed = 1;

//extern Int_t GlobalSeed;

static const Double_t PI=3.14159265358979312;

//static const bool     ConventionOld = true;
static const bool     ConventionOld = false;

class DIGIT_B2 : public TNamed {
 public:
  DIGIT_B2();
  DIGIT_B2(char *name, char *title);
  DIGIT_B2(char *name, char *title,
	   Int_t myglobalseed,
	   TString TransporModel = TString("HighResistivity"));
  
  DIGIT_B2(char *name, char *title,
	   Int_t myglobalseed,
	   Int_t aMapping,
	   Double_t aPitchU,   Double_t aPitchV,
	   Double_t ashiftU,   Double_t ashiftV,
	   Int_t    aNpixelsU, Int_t    aNpixelsV,
	   Double_t aEffectiveEpi,
	   Double_t aNoiseElectrons,
	   TString TransportModel);
  
  virtual ~DIGIT_B2();
 
  void         Run();
  Double_t GaussianLaw(Double_t mean, Double_t sigma); // OOO
  Int_t PoissonLaw(Float_t Lambda) ;
  Double_t LandauLaw(Double_t mean, Double_t sigma);

 
  Int_t GetPixelNumber(Double_t Xpos, Double_t Ypos);
  
  void GetXYPixelNumber(Int_t &Xpix, Int_t &Ypix,
			Int_t PixelNumber);
  
  void GetXYPixelCenter(Double_t &Xpix, Double_t &Ypix,
			Int_t PixelNumber);
  
  void Compute_CollectedCharge(Double_t SegmentSize,
			       Double_t InputUpos,Double_t InputVpos,
			       Double_t OutputUpos,Double_t OutputVpos,
			       Double_t Edeposited,
			       Double_t RandTheta,
			       Double_t RandPhi,
			       vector< Int_t > &fPixelMap,
			       vector< Float_t > &fAnalogChargeMap,
			       vector< Int_t > &fDigitalChargeMap);
  
  void Compute_digitisation(vector< Int_t > fPixelMap,
			    vector< Float_t > fAnalogChargeMap,
			    vector< Int_t > &fDigitalChargeMap,
			    Double_t Discri_Threshold_electrons);
  
  void  Transport_Noise_Digitisation(Double_t SegmentSize,
				     Double_t InputUpos,Double_t InputVpos,
				     Double_t OutputUpos,Double_t OutputVpos,
				     Double_t Edeposited,
				     Double_t RandTheta,
				     Double_t RandPhi,
				     vector< Int_t > &fPixelMap,vector< Float_t > &fAnalogChargeMap,
				     vector< Int_t > &fDigitalChargeMap,
				     Double_t Discri_Threshold_electrons);
  
  void      SetSensorParameters(Double_t  aPitchU,       Double_t  aPitchV,
				Double_t  ashiftU,       Double_t  ashiftV,
				Int_t     aNpixelsU,     Double_t  aNpixelsV,
				Double_t  aEffectiveEpi, Double_t  aNoiseElectrons,
				Int_t     aMapping,
				TString   aTransportModel);  //Sets the sensor parameters needed for the charge transport and digitization
  
  void      DefineTransportModel(void);  //Defines the Charge transport model using mainly the TheTransportModel model variable
  
  void      SetGlobalSeed(Long64_t  aSeed) {GlobalSeed = aSeed;}
  Long64_t  GetGlobalSeed()                {return  GlobalSeed;}
  
  Double_t  GetPitchU()                    {return  PitchU;}
  Double_t  GetPitchV()                    {return  PitchV;}
  Double_t  GetshiftU()                    {return  shiftU;}
  Double_t  GetshiftV()                    {return  shiftV;}
  Int_t     GetNpixelsU()                  {return  NpixelsU;}
  Int_t     GetNpixelsV()                  {return  NpixelsV;}
  Double_t  GetEffectiveEpi()              {return  EffectiveEpi;}
  Double_t  GetNoiseElectrons()            {return  NoiseElectrons;}
  Int_t     GetMapping()                   {return  Mapping;}
  TString   GetTransportModel()            {return  TheTransportModel;}
  bool      GetDoDigitization()            {return  DoDigitization;}
  
  //Same conventions as in TAF
  void      ComputePixelPositionUV_FromColRow(Int_t col, Int_t lin,   Double_t &u, Double_t &v);
  void      ComputePixelPositionColRow_FromUV(Double_t u, Double_t v, double &col, double &lin);
  
  void      ComputeUVRange(void);

 //protected: // OOO
  //ClassDef(DIGIT_B2,1)
 protected:
  Int_t Test;

  TF1 *mymodel1D_1st;
  TF1 *mymodel1D_2nd;
  TF2 *mymodel2D;

  Long64_t  GlobalSeed;
  Int_t     Mapping;
  Double_t  PitchU;
  Double_t  PitchV;
  Double_t  shiftU;
  Double_t  shiftV;
  Int_t     NpixelsU;
  Int_t     NpixelsV;
  Double_t  EffectiveEpi;
  Double_t  NoiseElectrons;
  TString   TheTransportModel;
  
  Double_t  MaxU;
  Double_t  MaxV;
  
  Double_t  RU[2];
  Double_t  RV[2];
  
  bool      DoDigitization;

};

//==============================================================================
#endif

