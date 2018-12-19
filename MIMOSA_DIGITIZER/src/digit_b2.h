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
#include "TStopwatch.h"
#include "TVector2.h"

using namespace std;

//Int_t GlobalSeed = 1;

//extern Int_t GlobalSeed;

static const Double_t PI = 3.14159265358979312;

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
	   Int_t aAnalysisMode,
	   Double_t aPitchU,   Double_t aPitchV,
	   Double_t ashiftU,   Double_t ashiftV,
	   Int_t    aNpixelsU, Int_t    aNpixelsV,
	   Double_t aEffectiveEpi,
	   Double_t aNoiseElectrons,
	   TString TransportModel,
	   Double_t aCalib       = 1.0,
	   Int_t    aADCbits     = 1,
	   Double_t aADCRangeMin = 0.0,
	   Double_t aADCRangeMax = 1.0);
  
  virtual ~DIGIT_B2();
 
  void         Run(Long64_t aGlobalSeed = 1);
  
  void         Run2(Long64_t    Nevents          = 1000,
                    Int_t       aMapping         = 1,
		    Int_t       aAnalysisMode    = 3,
		    Int_t       aNpixelsU        = 516,
		    Int_t       aNpixelsV        = 1152,
		    Double_t    aPitchU          = 18.4, //um
		    Double_t    aPitchV          = 18.4, //um
		    Double_t    aEffectiveEpi    = 10.0, //um
		    Double_t    aNoiseElectrons  = 14.0, //um
		    const char* aTransportModel  = "AMSlorgaus",
		    Double_t    Trk_thetaLoc     = 0.0,
		    Long64_t    aGlobalSeed      = 72984,
		    const char* output_file      = "output_file",
		    bool        verbose          = false);
  
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
				Int_t     aAnalysisMode,
				TString   aTransportModel,
				Double_t  aCalib       = 1.0,
				Int_t     aADCbits     = 1,
				Double_t  aADCRangeMin = 0.0,
				Double_t  aADCRangeMax = 1.0);  //Sets the sensor parameters needed for the charge transport and digitization
  
  void      DefineTransportModel(void);  //Defines the Charge transport model using mainly the TheTransportModel model variable
  
  void      SetGlobalSeed(Long64_t  aSeed) {GlobalSeed = aSeed;}
  Long64_t  GetGlobalSeed()                {return  GlobalSeed;}
  
  Double_t  GetPitchU()                    {return  PitchU;}
  Double_t  GetPitchV()                    {return  PitchV;}
  Double_t  GetshiftU()                    {return  shiftU;}
  Double_t  GetshiftV()                    {return  shiftV;}
  Int_t     GetNpixelsU()                  {return  NpixelsU;}
  Int_t     GetNpixelsV()                  {return  NpixelsV;}
  void      SetEffectiveEpi(Double_t aEpi) {EffectiveEpi = aEpi;}
  Double_t  GetEffectiveEpi()              {return  EffectiveEpi;}
  Double_t  GetNoiseElectrons()            {return  NoiseElectrons;}
  Int_t     GetMapping()                   {return  Mapping;}
  Int_t     GetAnalysisMode()              {return  AnalysisMode;}
  Double_t  GetCalib()                     {return  Calib;}
  Int_t     GetADCbits()                   {return  ADCbits;}
  Double_t  GetADCRangeMin()               {return  ADCRangeMin;}
  Double_t  GetADCRangeMax()               {return  ADCRangeMax;}
  TString   GetTransportModel()            {return  TheTransportModel;}
  bool      GetDoDigitization()            {return  DoDigitization;}
  
  //Same conventions as in TAF
  int       GetGlobalIdxFromColRow(int col, int row);
  void      GetColRowFromUV(Double_t u,Double_t v, int &col, int &row);
  void      ComputePixelPositionUV_FromColRow(Int_t col, Int_t lin,   Double_t &u, Double_t &v);
  void      ComputePixelPositionColRow_FromUV(Double_t u, Double_t v, double &col, double &lin);
  
  void      ComputeUVRange(void);
  
  int       RoundOff(double a);
  
  void      CheckADCRanges(void);
  
  void      Initialize(void);
  
  Double_t  ConvertToADC(Double_t Signal);
  
  double    GetTrkDistToDiode(double TrkU, double TrkV);
  TVector2  GetTrkShiftedPosition(double TrkU, double TrkV);

 //protected: // OOO
  //ClassDef(DIGIT_B2,1)
 protected:
  Int_t Test;

  TF1 *mymodel1D_1st;
  TF1 *mymodel1D_2nd;
  TF2 *mymodel2D;

  Long64_t  GlobalSeed;
  Int_t     Mapping;
  Int_t     AnalysisMode;
  Double_t  PitchU;
  Double_t  PitchV;
  Double_t  shiftU;
  Double_t  shiftV;
  Int_t     NpixelsU;
  Int_t     NpixelsV;
  Double_t  EffectiveEpi;
  Double_t  NoiseElectrons;
  Double_t  Calib;
  Int_t     ADCbits;
  Double_t  ADCRangeMin;
  Double_t  ADCRangeMax;
  TString   TheTransportModel;
  
  Double_t  MaxU;
  Double_t  MaxV;
  
  Double_t  RU[2];
  Double_t  RV[2];
  
  bool      DoDigitization;
  
  TStopwatch fWatch;
  Long64_t   fPrintFreq;

};

//==============================================================================
#endif

