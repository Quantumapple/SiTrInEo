#include <src/digit_b2.h>
#include <TStopwatch.h>
#include <TROOT.h> // for gROOT object
#include <TMath.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMarker.h>
#include <TF2.h>
#include <TF1.h>
#include <TVector2.h>
#include <TEllipse.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLegend.h>

#include <assert.h>

using namespace std;


float epsilon = 1.0e-3;

//==============================================================================
//ClassImp(DIGIT_B2)
DIGIT_B2::DIGIT_B2() : TNamed("DIGIT_B2","DIGIT_B2 title"), 
                       GlobalSeed(1),
                       Mapping(1),
                       AnalysisMode(3),
                       PitchU(18.4),  PitchV(18.4),
                       shiftU(0.0),   shiftV(0.0),
                       NpixelsU(576), NpixelsV(1152),
                       EffectiveEpi(10.0),
                       NoiseElectrons(14.0),
                       Calib(1.0),
                       ADCbits(1),
                       ADCRangeMin(0.0),
                       ADCRangeMax(1.0),
                       TheTransportModel("HR2016")
{
  //
  // default constructor
  Initialize();
  
  //
}
//______________________________________________________________________________
//  
DIGIT_B2::DIGIT_B2(char *name, char *title )  : 
                   TNamed(name,title),
                   GlobalSeed(1),
                   Mapping(1),
                   AnalysisMode(3),
                   PitchU(18.4),  PitchV(18.4),
                   shiftU(0.0),   shiftV(0.0),
                   NpixelsU(576), NpixelsV(1152),
                   EffectiveEpi(10.0),
                   NoiseElectrons(14.0),
                   Calib(1.0),
                   ADCbits(1),
                   ADCRangeMin(0.0),
                   ADCRangeMax(1.0),
                   TheTransportModel("HR2016")
{

  std::cout<<"+++++++++++++++++ START         +++++++++++++++++++++++++++++++"<<endl;
  Initialize();
  
  //Run();
  std::cout<<"+++++++++++++++++ END           +++++++++++++++++++++++++++++++"<<endl;

}
//______________________________________________________________________________
//  
DIGIT_B2::DIGIT_B2(char *name, char *title, Int_t myglobalseed, TString TransportModel)  : 
                   TNamed(name,title),
                   GlobalSeed(myglobalseed),
                   Mapping(1),
                   AnalysisMode(3),
                   PitchU(18.4),  PitchV(18.4),
                   shiftU(0.0),   shiftV(0.0),
                   NpixelsU(576), NpixelsV(1152),
                   EffectiveEpi(10.0),
                   NoiseElectrons(14.0),
                   Calib(1.0),
                   ADCbits(1),
                   ADCRangeMin(0.0),
                   ADCRangeMax(1.0),
                   TheTransportModel(TransportModel)
{

  std::cout<<"+++++++++++++++++ START         +++++++++++++++++++++++++++++++"<<endl;
  Initialize();
  
  // Run();
  std::cout<<"+++++++++++++++++ END           +++++++++++++++++++++++++++++++"<<endl;

}
//______________________________________________________________________________
//
DIGIT_B2::DIGIT_B2(char *name, char *title,
		   Int_t myglobalseed, 
		   Int_t aMapping,
		   Int_t aAnalysisMode,
		   Double_t aPitchU,   Double_t aPitchV,
		   Double_t ashiftU,   Double_t ashiftV,
		   Int_t    aNpixelsU, Int_t    aNpixelsV,
		   Double_t aEffectiveEpi,
		   Double_t aNoiseElectrons,
		   TString  TransportModel,
		   Double_t aCalib,
		   Int_t    aADCbits,
		   Double_t aADCRangeMin,
		   Double_t aADCRangeMax)  : 
		   TNamed(name,title),
		   GlobalSeed(myglobalseed),
		   Mapping(aMapping),
		   AnalysisMode(aAnalysisMode),
		   PitchU(aPitchU),     PitchV(aPitchV),
		   shiftU(ashiftU),     shiftV(ashiftV),
		   NpixelsU(aNpixelsU), NpixelsV(aNpixelsV),
		   EffectiveEpi(aEffectiveEpi),
		   NoiseElectrons(aNoiseElectrons),
		   Calib(aCalib),
		   ADCbits(aADCbits),
                   ADCRangeMin(aADCRangeMin),
                   ADCRangeMax(aADCRangeMax),
		   TheTransportModel(TransportModel)
{

  std::cout<<"+++++++++++++++++ START         +++++++++++++++++++++++++++++++"<<endl;
  
  Initialize();
  
  // Run();
  std::cout<<"+++++++++++++++++ END           +++++++++++++++++++++++++++++++"<<endl;

}
//______________________________________________________________________________
//
void  DIGIT_B2::SetSensorParameters(Double_t  aPitchU,       Double_t  aPitchV,
				    Double_t  ashiftU,       Double_t  ashiftV,
				    Int_t     aNpixelsU,     Double_t  aNpixelsV,
				    Double_t  aEffectiveEpi, Double_t aNoiseElectrons,
				    Int_t     aMapping,
				    Int_t     aAnalysisMode,
				    TString   aTransportModel,
				    Double_t  aCalib,
				    Int_t     aADCbits,
				    Double_t  aADCRangeMin,
				    Double_t  aADCRangeMax)
{
  
  Mapping           = aMapping;
  AnalysisMode      = aAnalysisMode;
  PitchU            = aPitchU;
  PitchV            = aPitchV;
  shiftU            = ashiftU;
  shiftV            = ashiftV;
  NpixelsU          = aNpixelsU;
  NpixelsV          = aNpixelsV;
  EffectiveEpi      = aEffectiveEpi;
  NoiseElectrons    = aNoiseElectrons;
  Calib             = aCalib;
  ADCbits           = aADCbits;
  ADCRangeMin       = aADCRangeMin;  // in Volts
  ADCRangeMax       = aADCRangeMax;  // in Volts
  TheTransportModel = aTransportModel;
  
  Initialize();
  
  return;
  
}
//______________________________________________________________________________
//
Int_t DIGIT_B2::PoissonLaw(Float_t Lambda) 
{ 
  Int_t n;
  GlobalSeed++;
  TRandom3 *r3 = new TRandom3(GlobalSeed);
  //Float_t smear = r3->Rndm(245);
  n = r3->Poisson(Lambda);

  //  TMath::PoissonI(x,Lambda)
  delete r3;
  return n;
    
}
//______________________________________________________________________________
//  
Double_t DIGIT_B2::LandauLaw(Double_t mean, Double_t sigma) 
{ 
  Double_t x;
  GlobalSeed++;
  TRandom3 *r3 = new TRandom3(GlobalSeed);
  //Float_t smear = r3->Rndm(245);
  x = r3->Landau(mean,sigma);
  //  TMath::PoissonI(x,Lambda)
  delete r3;
  return x;
}
//______________________________________________________________________________
//  
Double_t DIGIT_B2::GaussianLaw(Double_t mean, Double_t sigma) // XXXX
{ 
  Double_t x;
  GlobalSeed++;
  TRandom3 *r3 = new TRandom3(GlobalSeed);
  //Float_t smear = r3->Rndm(245);
  x = r3->Gaus(mean,sigma);
  //  TMath::PoissonI(x,Lambda)
  delete r3;
  return x;
}
//_______________________________________________________________________________________
//
  Double_t Lorentz2D(Double_t *x, Double_t *par){ 
    //x[0] = x
    //x[1] = y
    // par[0] = Gamma
    // par[1] = x0
    // par[2] = y0
    // par[3] = norm
  if(par[0]>0){
    Double_t Pi = 3.141592653;
    return par[3]*par[0]/Pi/((x[0]-par[1])*(x[0]-par[1])+(x[1]-par[2])*(x[1]-par[2])+par[0]*par[0]) ; 
  }else{
    return 0;
  }
}
//_______________________________________________________________________________________
//
Double_t SumGausLorentz2D(Double_t *x, Double_t *par){
  //par[0] Norm_1
  //par[1] x0_1
  //par[2] sigma_x_1
  //par[3] y0_1
  //par[4] sigma_y_1
  // par[5] = Gamma
  // par[6] = x0
  // par[7] = y0
  // par[8] = norm
  Double_t Pi = 3.141592653;
  if((par[2]!=0.0) && (par[4]!=0.0) ){
    double rx = (x[0]-par[1])/par[2];
    double ry = (x[1]-par[3])/par[4];
    return par[0]*( TMath::Exp(-(rx*rx+ry*ry)/2.0)
    		    +par[8]*par[5]/Pi/((x[0]-par[6])*(x[0]-par[6])+(x[1]-par[7])*(x[1]-par[7])+par[5]*par[5])
    		    );
  }else{
    return 0;  
  }
}
//_______________________________________________________________________________________
//
Double_t SumGausLorentz1D(Double_t *x, Double_t *par){
  //par[0] Norm_g
  //par[1] x0_g
  //par[2] sigma_g
  // par[3] = Gamma_lor
  // par[4] = x0_lor
  // par[5] = norm_lor
  /*   for (int j=0;j<10; j++) {
       cout<<"par " <<j<<" "<<par[j]<<endl;
       }
  */
  Double_t Pi = 3.141592653;
  if((par[2]!=0.0) ){
    Double_t rx = (x[0]-par[1])/par[2];
    Double_t tempoutput;
    tempoutput= par[0]*( TMath::Exp(-(rx*rx)/2.0)
		    +par[5]*par[3]/Pi/ ((x[0]-par[4])*(x[0]-par[4]) +par[3]*par[3])
    		    );
    //  cout<<"SumGausLorentz1D " <<tempoutput<<endl;
    return tempoutput;
   
  }else{
    return 0;
  }
}
//______________________________________________________________________________
//   
void  DIGIT_B2::ComputePixelPositionUV_FromColRow(Int_t col, Int_t lin, Double_t &u, Double_t &v)
{
  
  // Compute the 3D position of the strip at column "col" and line "lin",
  //  set the values in the "u, v, w" variables.
  // The strip position depends on the mapping, look below.
  //
  // BEWARE: if you add a new case here, add it as weel in ComputePixelPositionColRow_FromUV
  //
  
  switch (Mapping) {  
    case 1:
    default:
    {
      // When pixels are organized on an orthogonal grid
      u = ((2*col - NpixelsU + 1) * PitchU)/2;
      v = ((2*lin - NpixelsV + 1) * PitchV)/2;
      break;
    } 
    case 2:
    {
      // When pixels are staggered from one column to the other
      u = ((2*col - NpixelsU + 1 ) * PitchU)/2 ;
      if ( col%2 == 0 ) v = (((lin - NpixelsV/2.0 + 1./2. ) * PitchV)) + 0.30 * PitchV;  //clm 2012.11.24
      else              v = (((lin - NpixelsV/2.0 + 1./2. ) * PitchV)) - 0.19 * PitchV;  //clm 2012.11.24 
      break;
    } 
    case 3:
    {
      //Clm Mapping For M32 L8_2
      if ( lin==7 && col == 10) u = ((2*col - NpixelsU + 1 ) * PitchU)/2 ;  
      else           u = -9999; //((2*col+2 - NpixelsU + 1 ) * PitchU)/2 ;  
      v = ((2*lin - NpixelsV + 1 ) * PitchV)/2 ;
      break;
    } 
    case 4:
    {
      // MonteCarlo Simulation. // LC 2014/01/10.
      u = ((2.*col - NpixelsU + 2. ) * PitchU)/2. ;  // LC 2015/02/17 : +2 commented                                                                           
      v = ((2.*lin - NpixelsV + 2. ) * PitchV)/2. ;  // LC 2015/02/17 : +2 commented 
      break;
    } 
    case 5:
    {
      //Mapping FSBB
      u = (2*col + 1 - NpixelsU) * PitchU/2.0;
      double fraction = 0.25;
      //double fraction = 0.75;
      if ( col%2 == 0 ) v = 0.5*(NpixelsV - 2*lin - 2*(fraction    )) * PitchV;
      else              v = 0.5*(NpixelsV - 2*lin - 2*(1 - fraction)) * PitchV;
      break;
    } 
    case 6:
    {
      //Mapping FSBB-bis
      u = (2*col + 1 - NpixelsU) * PitchU/2.0;
      //double fraction = 0.25;
      double fraction = 0.75;
      if ( col%2 == 0 ) v = 0.5*(NpixelsV - 2*lin - 2*(fraction    )) * PitchV;
      else              v = 0.5*(NpixelsV - 2*lin - 2*(1 - fraction)) * PitchV;
      break;
    } 
    case 7:
    {
      //M22-THRB7 and B6 -> set the matrix type  and Mapping accordingly in the config file !
      u = ((2*col - NpixelsU + 1 ) * PitchU)/2 ;       
      if ( col%2 == 0 ) v = (((lin - NpixelsV/2.0 + 1./2. ) * PitchV)) - 0.50 * PitchV;  //correct +0.5
      else              v = (((lin - NpixelsV/2.0 + 1./2. ) * PitchV)) - 0.00 * PitchV;  //correct -0.5
      break;
    }
      
  }
  
  return;
  
}
//______________________________________________________________________________
//
void  DIGIT_B2::ComputePixelPositionColRow_FromUV(Double_t u, Double_t v, double &col, double &lin)
{
  

  // Compute the 2D position of a strip in the variables u and v to the set of variables column and line
  // The strip position depends on the mapping, look below.
  //
  // BEWARE: if you add a new case here, add it as well in ComputePixelPositionUV_FromColRow
  
  switch (Mapping) {
    case 1:
    default:
    {
      // When pixels are organized on an orthogonal grid
      col = (u/PitchU) + ((NpixelsU - 1)/2.);
      lin = (v/PitchV) + ((NpixelsV - 1)/2.);
      break;
    }
    case 2:
    {
      // When pixels are staggered from one column to the other
      col = (u/PitchU) + ((NpixelsU - 1)/2.);
      if(int(col)%2 == 0) lin = (v/PitchV) + ((NpixelsV - 1)/2.) - 0.30;
      else                lin = (v/PitchV) + ((NpixelsV - 1)/2.) + 0.19;
      break;
    }
    case 3:
    {
      //Clm Mapping For M32 L8_2
      col = (u/PitchU) + ((NpixelsU - 1)/2.);
      lin = (v/PitchV) + ((NpixelsV - 1)/2.);
      break;
    }
    case 4:
    {
      // MonteCarlo Simulation. // LC 2014/01/10.
      col = (u/PitchU) + ((NpixelsU /*- 2*/)/2.); // LC 2015/02/17 : -2 commented
      lin = (v/PitchV) + ((NpixelsV /*- 2*/)/2.); // LC 2015/02/17 : -2 commented
      break;
    } 
    case 5:
    {
      //Mapping FSBB
      col = (u/PitchU) + ((NpixelsU - 1)/2.);
      
      double fraction = 0.25;
      //double fraction = 0.75;
      
      int IntCol = int(col);
      double delta = col - IntCol;
      if(delta >= 0.5) IntCol++;
      
      if(IntCol%2 == 0) lin = ((NpixelsV - 2*(fraction    ))/2.0)  -  (v/PitchV);
      else              lin = ((NpixelsV - 2*(1 - fraction))/2.0)  -  (v/PitchV);
      
      break;
    }
    case 6:
    {
      //Mapping FSBB-bis
      col = (u/PitchU) + ((NpixelsU - 1)/2.);
      
      //double fraction = 0.25;
      double fraction = 0.75;
      
      int IntCol = int(col);
      double delta = col - IntCol;
      if(delta >= 0.5) IntCol++;
      
      if(IntCol%2 == 0) lin = ((NpixelsV - 2*(fraction    ))/2.0)  -  (v/PitchV);
      else              lin = ((NpixelsV - 2*(1 - fraction))/2.0)  -  (v/PitchV);
      
      break;
    } 
    case 7:
    {
      //M22-THRB7 and B6 -> set the matrix type and Mapping accordingly in the config file !
      col = (u/PitchU) + ((NpixelsU - 1)/2.);
      if(int(col)%2 == 0) lin = (v/PitchV) + ((NpixelsV - 1)/2.) + 0.50; //correct -0.5
      else                lin = (v/PitchV) + ((NpixelsV - 1)/2.) + 0.0; //correct +0.5
      break;
    }
      
  }
  
  return;
  
}
//______________________________________________________________________________
//   
// return the pixel number (integer) from the X and Y position (float)
Int_t DIGIT_B2::GetPixelNumber(Double_t Xpos,   Double_t Ypos)
{

  Int_t PixelNumber;
  
#if 0
  if(ConventionOld) {
    Int_t XPixelNumber = 0;
    Int_t YPixelNumber = 0;
  
    //compute the column odd or even:
    if( (int(Xpos/(PitchU)))%2 ==1){
      //odd column = no staggering
      //XPixelNumber = int(Xpos/(PitchU));
      //YPixelNumber = int(Ypos/(PitchV));
      XPixelNumber = RoundOff((Xpos/PitchU) - 0.5);
      YPixelNumber = RoundOff((Ypos/PitchV) - 0.5);
    }
    else{
      //even column = staggering
      //XPixelNumber = int((Xpos-shiftU) / PitchU );
      //YPixelNumber = int((Ypos-shiftV) / PitchV );
      XPixelNumber = RoundOff(((Xpos-shiftU)/PitchU) - 0.5);
      YPixelNumber = RoundOff(((Ypos-shiftV)/PitchV) - 0.5);
    }
    PixelNumber = XPixelNumber + NpixelsU*YPixelNumber;
  }
  else {
#endif
    double col,row;
    ComputePixelPositionColRow_FromUV(Xpos,Ypos,col,row);
    //PixelNumber = int(col) + NpixelsU*int(lin);
    PixelNumber = RoundOff(col) + NpixelsU*RoundOff(row);
  //}

  if(
    (Xpos  < RU[0] - epsilon)  ||
    (Xpos  > RU[1] + epsilon)  ||
    (Ypos  < RV[0] - epsilon)  ||
    (Ypos  > RV[1] + epsilon)
  ){
    cout<<" WARNING  charge is going outside the plane limits"<<endl;
    return 0;
  }
  else return PixelNumber;
  
}
//______________________________________________________________________________
//
void  DIGIT_B2::GetColRowFromUV(Double_t u,Double_t v, int &col, int &row)
{
  
  double the_col, the_row;
  ComputePixelPositionColRow_FromUV(u,v,the_col,the_row);
  col = RoundOff(the_col);
  row = RoundOff(the_row);
  
  return;
  
}
//______________________________________________________________________________
//   
// get Xnum(int) and Ynum(int) of the pixel from the pixel number (int)
void DIGIT_B2::GetXYPixelNumber(Int_t &Xpix, Int_t &Ypix,
				Int_t PixelNumber){

  Xpix = PixelNumber%(NpixelsU);
  Ypix = PixelNumber/(NpixelsU);
}
//______________________________________________________________________________
//   
// get the X and Y position (float) of the pixel from the pixel number (int), taking into account diode shift.
void DIGIT_B2::GetXYPixelCenter(Double_t &Xpix,  Double_t &Ypix,
				Int_t PixelNumber)
{
  
#if 0
  if(ConventionOld) {
    Xpix = (PitchV) * (0.5 + PixelNumber%NpixelsU)   + (PixelNumber%(NpixelsU)%2)*shiftU;
    Ypix = (PitchU) * (0.5 + PixelNumber/(NpixelsU)) + (PixelNumber%(NpixelsU)%2)*shiftV;
  }
  else {
#endif
    int col,lin;
    GetXYPixelNumber(col,lin,PixelNumber);
    ComputePixelPositionUV_FromColRow(col,lin,Xpix,Ypix);
  //}
  
  return;
  
}
//______________________________________________________________________________
//
int DIGIT_B2::GetGlobalIdxFromColRow(int col, int row)
{
  
  return col + NpixelsU*row;
  
}
//______________________________________________________________________________
//
void  DIGIT_B2::DefineTransportModel(void)
{
  
  DoDigitization = true;
  
  mymodel1D_1st = NULL;
  mymodel1D_2nd = NULL;
  
  //-----------------------------------------------------
  //-------------------------Compute charge transport up to the diodes
  mymodel1D_1st = NULL;
  mymodel1D_2nd = NULL;
  mymodel2D     = NULL;
  //-----------------------------------------------------
  //--ChargeModel 3
  //-----------------------------------------------------
  //---Lor+gaus Model
  Double_t LorGaussModel_Norm1_Cp0 =  0.160573;
  Double_t LorGaussModel_Norm1_Cp1 = -0.00184473;
  Double_t LorGaussModel_Norm1_Cp2 =  5.07964e-05;
  Double_t LorGaussModel_sigma_Cp0 =  0.95;
  Double_t LorGaussModel_sigma_Cp1 =  0.60;
  Double_t LorGaussModel_C_Cp0     =  0.171697;
  Double_t LorGaussModel_C_Cp1     =  0.316165;
  Double_t LorGaussModel_Norm_Cp0  =  5.69809;
  Double_t LorGaussModel_Norm_Cp1  =  2.55373;


  //-----------------------------------------------------
  //--ChargeModel 5
  //-----------------------------------------------------
  // 1 dimension model with 2 functions for each squares around the hit
  //-----------------------------------------------------
  Double_t Norm_g_1st     =  0.458955;
  Double_t x0_g_1st       = -3.98149;
  Double_t sigma_g_1st    = 13.1559;
  Double_t Gamma_lor_1st  =  3.98673;
  Double_t x0_lor_1st     =  1.79712;
  Double_t norm_lor_1st   =  6.4533;
  
  Double_t Norm_g_2nd     =  0.116703;
  Double_t x0_g_2nd       = -1.0688;
  Double_t sigma_g_2nd    = 17.4823;     
  Double_t Gamma_lor_2nd  = 47.0837;
  Double_t x0_lor_2nd     = -4.63879;
  Double_t norm_lor_2nd   =  3.71411;
  //-----------------------------------------------------
  
 //-----------------------------------------------------
  //--ChargeModel 6
  //-----------------------------------------------------
  Double_t A_sigma_g_1st     = 0.025;
  Double_t B_sigma_g_1st     = 0.5;
  Double_t A_Gamma_lor_1st   = 0.025;
  Double_t B_Gamma_lor_1st   = 0.5  ;
  Double_t A_sigma_g_2nd     = 0.025;
  Double_t B_sigma_g_2nd     = 0.5;
  Double_t A_Gamma_lor_2nd   = 0.025;
  Double_t B_Gamma_lor_2nd   = 0.5 ;
    
  Double_t sigma_1st = sigma_g_1st   * (A_sigma_g_1st   * PitchU + B_sigma_g_1st);
  Double_t Gamma_1st = Gamma_lor_1st * (A_Gamma_lor_1st * PitchU + B_Gamma_lor_1st);
  Double_t sigma_2nd = sigma_g_2nd   * (A_sigma_g_2nd   * PitchU + B_sigma_g_2nd);
  Double_t Gamma_2nd = Gamma_lor_2nd * (A_Gamma_lor_2nd * PitchU + B_Gamma_lor_2nd);

  //-----------------------------------------------------
  // CHOSE MODEL:
  //-----------------------------------------------------
  if((TheTransportModel == TString("AMSlorgaus")) || (TheTransportModel == TString("3"))) {
    //-----------------------------------------------------
    // Lorentz + Gauss model 
    //-----------------------------------------------------
    //par[0] = Norm_1
    //par[1] = x0_1
    //par[2] = sigma_x_1
    //par[3] = y0_1
    //par[4] = sigma_y_1
    //par[5] = Gamma
    //par[6] = x0
    //par[7] = y0
    //par[8] = norm    
    Double_t Norm1Value = LorGaussModel_Norm1_Cp0 + PitchU * LorGaussModel_Norm1_Cp1 + PitchU*PitchU*LorGaussModel_Norm1_Cp2;
    Double_t sigmaValue = LorGaussModel_sigma_Cp0 + PitchU * LorGaussModel_sigma_Cp1;
    Double_t Cvalue     = LorGaussModel_C_Cp0     + PitchU * LorGaussModel_C_Cp1;
    Double_t NormValue  = LorGaussModel_Norm_Cp0  + PitchU * LorGaussModel_Norm_Cp1;
    Double_t rangelimit_inpitchunit = 2.5;
    if(!mymodel2D) {
      mymodel2D = new TF2("funlor2d",SumGausLorentz2D,
			  -rangelimit_inpitchunit*PitchU,rangelimit_inpitchunit*PitchU,
			  -rangelimit_inpitchunit*PitchV,rangelimit_inpitchunit*PitchV,
			  9);
      mymodel2D->SetParNames("Norm_1","x0_1","sigma_x_1","y0_1","sigma_y_1","Gamma","x0","y0","norm");
    }
    else {
      mymodel2D->SetRange(-rangelimit_inpitchunit*PitchU,rangelimit_inpitchunit*PitchU,
			  -rangelimit_inpitchunit*PitchV,rangelimit_inpitchunit*PitchV);
    }
    Double_t params1[] = {Norm1Value,0.0,sigmaValue,0.0,sigmaValue,Cvalue,0.0,0.0,NormValue};
    mymodel2D->SetParameters(params1);
  }
  //-----------------------------------------------------
  else if((TheTransportModel == TString("HR2015")) || (TheTransportModel == TString("5"))) {
    //--RangeLimit_InPitchUnit 2.5
    if(!mymodel1D_1st) {
      mymodel1D_1st = new TF1("namemodel1D_1st",SumGausLorentz1D,0,4.0*PitchU,6);
      mymodel1D_1st->SetParNames("Norm_g","x0_g","sigma_g","Gamma_lor","x0_lor","norm_lor");
    }
    else {
      mymodel1D_1st->SetRange(0,4.0*PitchU);
    }
    
    if(!mymodel1D_2nd) {
      mymodel1D_2nd = new TF1("namemodel1D_2nd",SumGausLorentz1D,0,4.0*PitchU,6);
      mymodel1D_2nd->SetParNames("Norm_g","x0_g","sigma_g","Gamma_lor","x0_lor","norm_lor");
    }
    else {
      mymodel1D_2nd->SetRange(0,4.0*PitchU);
    }
    
    Double_t params1[] = {Norm_g_1st,x0_g_1st,sigma_g_1st, Gamma_lor_1st,x0_lor_1st,norm_lor_1st};
    Double_t params2[] = {Norm_g_2nd,x0_g_2nd,sigma_g_2nd, Gamma_lor_2nd,x0_lor_2nd,norm_lor_2nd};
    mymodel1D_1st->SetParameters(params1);
    mymodel1D_2nd->SetParameters(params2);
  }
  //-----------------------------------------------------
  else if((TheTransportModel == TString("HR2016")) || (TheTransportModel == TString("6"))) {
    if(!mymodel1D_1st) {
      mymodel1D_1st = new TF1("namemodel1D_1st",SumGausLorentz1D,0,4.0*PitchU,6);
      mymodel1D_1st->SetParNames("Norm_g","x0_g","sigma_g","Gamma_lor","x0_lor","norm_lor");
    }
    else {
      mymodel1D_1st->SetRange(0,4.0*PitchU);
    }
    
    if(!mymodel1D_2nd) {
      mymodel1D_2nd = new TF1("namemodel1D_2nd",SumGausLorentz1D,0,4.0*PitchU,6);
      mymodel1D_2nd->SetParNames("Norm_g","x0_g","sigma_g","Gamma_lor","x0_lor","norm_lor");
    }
    else {
      mymodel1D_2nd->SetRange(0,4.0*PitchU);
    }
    
    Double_t params1[] = {Norm_g_1st,x0_g_1st,sigma_1st, Gamma_1st,x0_lor_1st,norm_lor_1st};
    Double_t params2[] = {Norm_g_2nd,x0_g_2nd,sigma_2nd, Gamma_2nd,x0_lor_2nd,norm_lor_2nd};
    mymodel1D_1st->SetParameters(params1);
    mymodel1D_2nd->SetParameters(params2);
  }
  else { //default: model 5
    if(!mymodel1D_1st) {
      mymodel1D_1st = new TF1("namemodel1D_1st",SumGausLorentz1D,0,4.0*PitchU,6);
      mymodel1D_1st->SetParNames("Norm_g","x0_g","sigma_g","Gamma_lor","x0_lor","norm_lor");
    }
    else {
      mymodel1D_1st->SetRange(0,4.0*PitchU);
    }
    
    if(!mymodel1D_2nd) {
      mymodel1D_2nd = new TF1("namemodel1D_2nd",SumGausLorentz1D,0,4.0*PitchU,6);
      mymodel1D_2nd->SetParNames("Norm_g","x0_g","sigma_g","Gamma_lor","x0_lor","norm_lor");
    }
    else {
      mymodel1D_2nd->SetRange(0,4.0*PitchU);
    }
    
    Double_t params1[] = {Norm_g_1st,x0_g_1st,sigma_g_1st, Gamma_lor_1st,x0_lor_1st,norm_lor_1st};
    Double_t params2[] = {Norm_g_2nd,x0_g_2nd,sigma_g_2nd, Gamma_lor_2nd,x0_lor_2nd,norm_lor_2nd};
    mymodel1D_1st->SetParameters(params1);
    mymodel1D_2nd->SetParameters(params2);
    
    cout << endl;
    cout << "Digitization Model " << TheTransportModel.Data() << "  doesn't correspond to list:" << endl;
    cout << " - AMSlorgaus or 3: Some explanation" << endl;
    cout << " - HR2015     or 5: Some explanation" << endl;
    cout << " - HR2016     or 6: Some explanation" << endl;
    cout << endl;

    DoDigitization = false;
  }
  
  return;
  
}
//______________________________________________________________________________
//
void DIGIT_B2::Compute_CollectedCharge(Double_t SegmentSize,
				       Double_t InputUpos,Double_t InputVpos,
				       Double_t OutputUpos,Double_t OutputVpos,
				       Double_t Edeposited,
				       Double_t RandTheta,
				       Double_t RandPhi,
				       vector< Int_t > &fPixelMap,vector< Float_t > &fAnalogChargeMap, 
				       vector< Int_t > &fDigitalChargeMap){
  //start
  //cout<<" DEBUG 100"<<endl;
  
  /*cout << "NpixelsU = " << NpixelsU << "  NpixelsV = " << NpixelsV << endl;
  cout << "PitchU = " << PitchU << "  PitchV = " << PitchV << endl;
  cout << "MaxU = " << MaxU << "  MaxV = " << MaxV << endl;*/
  
  //-----------------------------------------------------
  //--------------------------Compute energy deposition
  //-----------------------------------------------------
  //Double_t totalUlength=EffectiveEpi *TMath::Tan(RandTheta)*TMath::Cos(RandPhi);
  //Double_t totalVlength=EffectiveEpi *TMath::Tan(RandTheta)*TMath::Sin(RandPhi);
  Double_t totallentgh = sqrt(pow(InputUpos - OutputUpos,2) + 
                              pow(InputVpos - OutputVpos,2) + 
                              pow(EffectiveEpi,2));
  Float_t EnergyMPV   = 800.0 * totallentgh / 10.0;   
  Float_t EnergySIGMA = 180.0 * totallentgh / 10.0;

  //Used input deposited energy by default
  Float_t Energy = Edeposited;
  if(Energy < 0) {
    //If input deposited energy is negative generate deposited energy with Landau distribution
    Energy = LandauLaw(EnergyMPV,EnergySIGMA); //deposited energy in electrons.
    
    while(Energy > 20000) {
      Energy = LandauLaw(EnergyMPV,EnergySIGMA);
    }
  }
  //-----------------------------------------------------
  //-------------------------Compute the segment and energy on each segment
  // simplified model: here the segments have a fixed size.
  //-----------------------------------------------------
 
  Int_t fNSegment; //number of segments
  Double_t ChargePerSegment; //charge in each segement (in e-)
  if((Energy > 0) && (SegmentSize > 0.0)){
    fNSegment        = int(totallentgh*1.000001/SegmentSize) ;  //number of segments
    ChargePerSegment = Energy / float(fNSegment);
  }
  else {
    Energy           = 0.0;
    ChargePerSegment = 0;
    fNSegment        = 1;
  }
  vector<Float_t> fSegmentX;
  vector<Float_t> fSegmentY;
  vector<Float_t> fSegmentZ;
  vector<Float_t> fSegmentCharge;

  Double_t xstep = OutputUpos - InputUpos;
  Double_t ystep = OutputVpos - InputVpos;
  Double_t zstep = EffectiveEpi;
  //  cout<<" DEBUG 300"<<endl;
  
  for (Int_t i=0 ; i<fNSegment ; i++){
    fSegmentX.push_back(InputUpos + (float(i+0.5)* xstep/float(fNSegment)) );
    fSegmentY.push_back(InputVpos + (float(i+0.5)* ystep/float(fNSegment)) );
    fSegmentZ.push_back(            (float(i+0.5)* zstep/float(fNSegment)) );
    fSegmentCharge.push_back(ChargePerSegment);
  }

  // cout<<" DEBUG 400"<<endl;
  //--------loop on segments
  for (Int_t i=0; i<fNSegment; i++){
    int col,row;
    GetColRowFromUV(fSegmentX[i],fSegmentY[i],col,row);
    
    Double_t  xdpos, ydpos;
    ComputePixelPositionUV_FromColRow(col,row,xdpos, ydpos);
    
    Double_t TotalProba = 0.0;
    const Int_t Npix(25);
    //Double_t PixelposX[Npix];
    //Double_t PixelposY[Npix];
    Double_t Pixelproba[Npix];
    Int_t    PixelNumber[Npix];
    for(Int_t j=0; j<Npix; j++){
      //PixelposX[j]   = 0.0;
      //PixelposY[j]   = 0.0;
      Pixelproba[j]  = 0.0;
      PixelNumber[j] = 0;
    }
    
    for (Int_t j=0; j<Npix; j++){
      Int_t delta_col = (j%5) - 2;
      Int_t delta_row = (j/5) - 2;
      Float_t xeval,yeval;
      
      Double_t  xdposj, ydposj;
      int colj = col + delta_col;
      int rowj = row + delta_row;
      ComputePixelPositionUV_FromColRow(colj,rowj,xdposj, ydposj);
      
      xeval  = -fSegmentX[i] + xdposj;
      yeval  = -fSegmentY[i] + ydposj;
      Float_t xydist = sqrt(pow(xeval,2) + pow(yeval,2));
      
      //find if pixel is one of the four in the first square around the hit:
      if((TheTransportModel == TString("1")) || (TheTransportModel == TString("AMSlorgaus")) || (TheTransportModel == TString("3"))) {
	Pixelproba[j] = mymodel2D->Eval(xeval,yeval);
      }
      else {
	if((fabs(xeval) <= PitchU) && (fabs(yeval) <= PitchV)) Pixelproba[j] = mymodel1D_1st->Eval(xydist);
	else                                                   Pixelproba[j] = mymodel1D_2nd->Eval(xydist);
      }

      //PixelposX[j] = xdposj;
      //PixelposY[j] = ydposj;
      
      if( (colj >= 0 && colj <= NpixelsU-1) && (rowj >= 0 && rowj <= NpixelsV-1) ) PixelNumber[j] = GetGlobalIdxFromColRow(colj,rowj);
      else                                                                         PixelNumber[j] = -1;
      
      //if( (PixelposX[j] < RU[0]) || (PixelposX[j] > RU[1]) || (PixelposY[j] < RV[0]) || (PixelposY[j] > RV[1]) ) {
      //  PixelNumber[j] = -1; //outside the matrix.
      //}
      //else PixelNumber[j] = GetPixelNumber(PixelposX[j],PixelposY[j]);
      
      TotalProba  += Pixelproba[j];
    }
    
    //---loop on the 25 pixels around the seed to add:
    //   1. the pixels indices to fPixelMap
    //   2. the charge to the fAnalogChargeMap
    //   3. nothing to the fDigitalChargeMap (just reserve the memory)
    
    for (Int_t j=0; j<Npix; j++){
      Double_t AnalogCharge = 0.0;
      if(PixelNumber[j] >= 0){
	AnalogCharge = fSegmentCharge[i]*Pixelproba[j]/TotalProba;
	Bool_t found = false;
	Int_t k=0; // j = indice de la matrice 5x5  / k = indice sur fPixelMap / i = indice sur les segments
	while((!found) && (k<int(fPixelMap.size()))){
	  if(PixelNumber[j] == fPixelMap[k]){  // the pixel already exists in the list.   
	    found = true;
	    fAnalogChargeMap[k] += AnalogCharge; //just update the analog charge.
	  }
	  k++;
	}
	
	if(!found){  // the pixel doesn't exists in the list, so add it.
	  //fNpixels++;
	  fPixelMap.push_back(PixelNumber[j]);
	  fAnalogChargeMap.push_back(AnalogCharge);
	  fDigitalChargeMap.push_back(0);
	}
	
      }
    } //--------END loop on Npix
  } //--------END loop on segments

  return;
  
}
//______________________________________________________________________________
//  
void DIGIT_B2::Compute_digitisation(vector< Int_t > fPixelMap,vector< Float_t > fAnalogChargeMap, vector< Int_t > &fDigitalChargeMap, 
				    Double_t Discri_Threshold_electrons ){

  for (Int_t i = 0; i < int(fPixelMap.size()) ; i++){
    if (fAnalogChargeMap[i] <= Discri_Threshold_electrons){
      fDigitalChargeMap[i]=0;
    }else{
      fDigitalChargeMap[i]=1;
    }
  }
}
//______________________________________________________________________________
//  
void  DIGIT_B2::Transport_Noise_Digitisation(Double_t SegmentSize,
					     Double_t InputUpos,Double_t InputVpos,
					     Double_t OutputUpos,Double_t OutputVpos,
					     Double_t Edeposited,
					     Double_t RandTheta,
					     Double_t RandPhi,
					     vector< Int_t > &fPixelMap,vector< Float_t > &fAnalogChargeMap,
					     vector< Int_t > &fDigitalChargeMap,
					     Double_t Discri_Threshold_electrons)

{

  //----------compute collected charge:
  Compute_CollectedCharge(SegmentSize,
			  InputUpos,InputVpos,
			  OutputUpos,OutputVpos,
			  Edeposited,
			  RandTheta,
			  RandPhi,
			  fPixelMap,fAnalogChargeMap,fDigitalChargeMap);
  //------very simple way to add random noise (just gaussian shape)
  Double_t Noise;
  for (Int_t i = 0; i < int(fPixelMap.size()) ; i++){
    Noise =   GaussianLaw(0.0, NoiseElectrons);
    fAnalogChargeMap[i]+=Noise;
  }
  
  //  ----now do the Analog to Digital conversion    Compute_digitisation();
  Compute_digitisation(fPixelMap,fAnalogChargeMap,fDigitalChargeMap,Discri_Threshold_electrons);
  
  return;

}

//______________________________________________________________________________
//
void DIGIT_B2::ComputeUVRange(void)
{
  
  MaxU  = NpixelsU*PitchU;
  MaxV  = NpixelsV*PitchV;
  
  RU[0] = RV[0] = +1.0e+20;
  RU[1] = RV[1] = -1.0e+20;
  
  int col,lin;
  double U,V;
  
  for(int icol=0;icol<2;icol++) {
    for(int ilin=0;ilin<2;ilin++) {
      col = icol;
      lin = ilin;
      ComputePixelPositionUV_FromColRow(col,lin,U,V);
      U -= PitchU*0.5;
      V -= PitchV*0.5;
      if(RU[0] > U) RU[0] = U;
      if(RU[1] < U) RU[1] = U;
      if(RV[0] > V) RV[0] = V;
      if(RV[1] < V) RV[1] = V;
      
      col = NpixelsU - 1 - icol;
      lin = NpixelsV - 1 - ilin;
      ComputePixelPositionUV_FromColRow(col,lin,U,V);
      U += PitchU*0.5;
      V += PitchV*0.5;
      if(RU[0] > U) RU[0] = U;
      if(RU[1] < U) RU[1] = U;
      if(RV[0] > V) RV[0] = V;
      if(RV[1] < V) RV[1] = V;
    }
  }
  
  return;
  
}
//______________________________________________________________________________
//
int  DIGIT_B2::RoundOff(double a)
{
  
  int a_int = int(a);
  double delta = a - a_int;
  if(delta < 0.5) return a_int;
  else            return a_int + 1;
  
}
//______________________________________________________________________________
//
void DIGIT_B2::CheckADCRanges(void)
{
  
  if(AnalysisMode == 2) {
    if(ADCRangeMin >= ADCRangeMax) {
      cout << endl;
      cout << " ERROR: ADCRangesMax <= ADCRangesMin. Check your inputs. Exiting now!!!" << endl;
      cout << endl;
      
      assert(false);
    }
  }
  
  return;
  
}
//______________________________________________________________________________
//
void DIGIT_B2::Initialize(void)
{
  
  CheckADCRanges();
  ComputeUVRange();
  DefineTransportModel();
  
  return;
  
}
//______________________________________________________________________________
//
Double_t DIGIT_B2::ConvertToADC(Double_t Signal)
{
 
  Double_t ADCSignal = Signal;
  
  return ADCSignal;
  
}
//______________________________________________________________________________
//
double    DIGIT_B2::GetTrkDistToDiode(double tu, double tv)
{ 

  double distance = 1.0e+20;
  int col,lin;
  double col_p,lin_p;
  double u,v;

  //Find the pixel which the track intersects
  ComputePixelPositionColRow_FromUV(tu,tv,col_p,lin_p);
  col = RoundOff(col_p);
  lin = RoundOff(lin_p);

  for(int iii=-1;iii<2;iii++) {
    //Looping on the colums to the left and to the right of the main pixel 
    int col_tmp = col + iii;

    //Cut to ensure that the columns tested are inbetween 0 - Ncolumns-1
    if(col_tmp < 0) continue;
    if(col_tmp > GetNpixelsU()-1) continue;
    for(int jjj=-1;jjj<2;jjj++) {
      int lin_tmp = lin + jjj;

      //Cut to ensure that the lines tested are inbetween 0 - Nlines-1
      if(lin_tmp < 0) continue;
      if(lin_tmp >  GetNpixelsV()-1) continue;
      
      ComputePixelPositionUV_FromColRow(col_tmp,lin_tmp,u,v);
      float distance_tmp = sqrt(pow(tu-u,2) + pow(tv-v,2));

      if(distance > distance_tmp) distance = distance_tmp;

    }
  }
  
  return distance;
  
}
//______________________________________________________________________________
//
TVector2  DIGIT_B2::GetTrkShiftedPosition(double tu, double tv)
{
  
  double u = (tu + 0.5 * GetNpixelsU() * GetPitchU())/(2.0 * GetPitchU());
  u = (u - int(u))*2.0 * GetPitchU();
  double v = (tv + 0.5 * GetNpixelsV() * GetPitchV())/(2.0 * GetPitchV());
  v = (v - int(v))*2.0*GetPitchV();
  
  TVector2 Trk_Shifted_Pos(u,v);
   
  return Trk_Shifted_Pos;
  
}
//______________________________________________________________________________
//  
DIGIT_B2::~DIGIT_B2() { // 
  // virtual destructor
  //
  
  delete  mymodel1D_1st;
  delete  mymodel1D_2nd;
  delete  mymodel2D;
  
}
//______________________________________________________________________________
//
void  DIGIT_B2::Run(Long64_t aGlobalSeed) 
{ 
  std::cout<<"DIGIT_B2::Run() ..."<<endl;

  //Define Pixel Matrix:
  // zero position is assumed to be on pixel (0,0).
  // Theta is the polar angle (theta = 0 -> orthogonal track ; tetha ~ pi/2 -> low incident angle track
  // phi is the azimuthal angle (phi = 0 ->  track in the U direction).
  // all dimensions are in microns
  
  SetGlobalSeed(aGlobalSeed);
  
  Int_t    aMapping        = 1;
  Int_t    aAnalysisMode   = 3;
  Int_t    aNpixelsU       = 576;
  Int_t    aNpixelsV       = 1152;
  Double_t aPitchU         = 18.4;
  Double_t aPitchV         = 18.4;
  Double_t ashiftU         = 0; // for staggered pixel diodes
  Double_t ashiftV         = 0; // for staggered pixel diodes
  Double_t aEffectiveEpi   = 10.0;
  Double_t aNoiseElectrons = 14.0; // noise in electrons 
  TString aTransportModel("AMSlorgaus");
  //TString aTransportModel("HR2016");

  SetSensorParameters(aPitchU,       aPitchV,
                      ashiftU,       ashiftV,
                      aNpixelsU,     aNpixelsV,
                      aEffectiveEpi, aNoiseElectrons,
		      aMapping,
		      aAnalysisMode,
		      aTransportModel);
  
  Double_t Discri_Threshold_Noiseunits = 6.0;
  Double_t Discri_Threshold_electrons = Discri_Threshold_Noiseunits*NoiseElectrons;


  Double_t ThetaMaxdeg = 88.0;
  Double_t ThetaMaxrad = ThetaMaxdeg*PI/180.0;

  Double_t SegmentSize = 1.0; // size of segment for charge deposition in the epi layer
  
  Int_t Nevents = 20;
  
  TRandom3 *r3 = new TRandom3(GlobalSeed);

  TH1F *hmultiplicity = new TH1F("hmultiplicity","cluster multiplicity ",40,0,40);

  //----------loop on events:

  for(Int_t Nev = 0; Nev <Nevents; Nev++){
    //Generate random hit:
    //GlobalSeed++;
    //r3->SetSeed(GlobalSeed);
    
    Double_t RandU     = r3->Rndm()*(RU[1] - RU[0]) + RU[0];
    Double_t RandV     = r3->Rndm()*(RV[1] - RV[0]) + RV[0];
    Double_t RandTheta = r3->Rndm()*ThetaMaxrad;//0.0;//r3->Rndm()*ThetaMaxrad;
    Double_t RandPhi   = r3->Rndm()*2.0*PI;
    
    Double_t InputUpos    = RandU;
    Double_t InputVpos    = RandV;
    Double_t totalUlength = EffectiveEpi*TMath::Tan(RandTheta)*TMath::Cos(RandPhi);
    Double_t totalVlength = EffectiveEpi*TMath::Tan(RandTheta)*TMath::Sin(RandPhi);
    Double_t OutputUpos   = InputUpos + totalUlength;
    Double_t OutputVpos   = InputVpos + totalVlength;
        
    //pixel with deposited charge list
    //Int_t fNpixels = 0;
    vector< Int_t >   fPixelMap;          // index of pixels which collects some charge
    vector< Float_t > fAnalogChargeMap;   // Analog chargd in e- units
    vector< Int_t >   fDigitalChargeMap;  // digital charge in ADC units
    
    Transport_Noise_Digitisation(SegmentSize,
				 InputUpos,InputVpos,
				 OutputUpos,OutputVpos,
				 -1.0,
				 RandTheta,
				 RandPhi,
				 fPixelMap,fAnalogChargeMap,fDigitalChargeMap,
				 Discri_Threshold_electrons);


    //---print the output of the event:
    Int_t totalMult=0;
    for(Int_t i = 0; i < int(fPixelMap.size()) ; i++){
      Int_t nU,nV;
      Double_t U,V;
      GetXYPixelNumber(nU,nV,fPixelMap[i]);
      GetXYPixelCenter(U,V,fPixelMap[i]);
      
      totalMult+=fDigitalChargeMap[i];
    }
    hmultiplicity->Fill(totalMult);

  } // --------------END LOOP ON EVENTS
  
  hmultiplicity->Draw();
  
  delete r3;
  
  return;
  
}
//______________________________________________________________________________
//
void  DIGIT_B2::Run2(Long64_t    Nevents,
                     Int_t       aMapping,
		     Int_t       aAnalysisMode,
		     Int_t       aNpixelsU,
		     Int_t       aNpixelsV,
		     Double_t    aPitchU,
		     Double_t    aPitchV,
		     Double_t    aEffectiveEpi,
		     Double_t    aNoiseElectrons,
		     const char* aTransportModel,
		     Double_t    Trk_thetaLoc,
		     Long64_t    aGlobalSeed,
		     const char* output_file,
		     bool verbose)
{ 
  
  //Define Pixel Matrix:
  // zero position is assumed to be on pixel (0,0).
  // Theta is the polar angle (theta = 0 -> orthogonal track ; tetha ~ pi/2 -> low incident angle track
  // phi is the azimuthal angle (phi = 0 ->  track in the U direction).
  // all dimensions are in microns

  SetGlobalSeed(aGlobalSeed);
  
  std::vector< Int_t >   fPixelMap;          // index of pixels which collects some charge
  std::vector< Float_t > fAnalogChargeMap;   // Analog chargd in e- units
  std::vector< Int_t >   fDigitalChargeMap;  // digital charge in ADC units    
  
  const Int_t Nhits_check(20);
  double R_Thr[2];
  R_Thr[0] =  3.0 - 0.5;
  R_Thr[1] = 12.0 + 0.5;
  const int N_Thr_steps(R_Thr[1] - R_Thr[0]);
  
  Double_t ashiftU = 0;
  Double_t ashiftV = 0;

  SetSensorParameters(aPitchU,       aPitchV,
                      ashiftU,       ashiftV,
                      aNpixelsU,     aNpixelsV,
                      aEffectiveEpi, aNoiseElectrons,
		      aMapping,
		      aAnalysisMode,
		      TString(aTransportModel));
  
  //Double_t  Discri_Threshold_electrons = Thre_Noise*NoiseElectrons;
  
  Double_t SegmentSize = GetEffectiveEpi()/20.0; // size of segment for charge deposition in the epi layer

  TRandom3 *r3 = new TRandom3(GlobalSeed);
  
  if(verbose) {
    cout << endl;
    cout << "Setting-up digitization parameters:" << endl;
    cout << " - PitchU          = " << GetPitchU()         << " um"        << endl;
    cout << " - PitchV          = " << GetPitchV()         << " um"        << endl;
    cout << " - NpixU           = " << GetNpixelsU()                       << endl;
    cout << " - NpixV           = " << GetNpixelsV()                       << endl;
    cout << " - Epu-layer       = " << GetEffectiveEpi()   << " um"        << endl;
    cout << " - Noise           = " << GetNoiseElectrons() << " electrons" << endl;
    cout << " - Mapping         = " << GetMapping()                        << endl;
    cout << " - AnalysisMode    = " << GetAnalysisMode()                   << endl;
    cout << " - transport model = " << GetTransportModel().Data()          << endl;
    cout << " - SegmentSize     = " << SegmentSize         << " um"        << endl;
    cout << " - Trk theta_Loc   = " << Trk_thetaLoc        << " deg"       << endl;
    cout << " - GlobalSeed      = " << GetGlobalSeed()                     << endl;
    
    //cout << " - Threshold       = " << Thre_Noise << "/" << Discri_Threshold_electrons << " xNoise/electrons" << endl;
    
    cout << endl;
  }
  
  fPrintFreq = 10000;
  
  double TheSize = 0.05;
  
  gStyle->SetOptFit(1112);
  TString FitModel("([2]/(sqrt(2.0*TMath::Pi())*[1]))*exp( -0.5*pow((x - [0])/[1],2))");
  
  const int MaxPixelInCluster(200);
  int   NhitsToDraw[N_Thr_steps];
  TH2F* h_ClusterShape_Analog[Nhits_check][N_Thr_steps];
  TH2F* h_ClusterShape_Digital[Nhits_check][N_Thr_steps];
  TEllipse* InPostion[Nhits_check][N_Thr_steps];
  TEllipse* OutPostion[Nhits_check][N_Thr_steps];
  TEllipse* InPostionUV[Nhits_check][N_Thr_steps];
  TEllipse* OutPostionUV[Nhits_check][N_Thr_steps];
  int       NPixelsCheck[Nhits_check][N_Thr_steps];
  TVector2  PixelPosition[Nhits_check][N_Thr_steps][MaxPixelInCluster];
  double    PixelAnalogCharge[Nhits_check][N_Thr_steps][MaxPixelInCluster];
  double    PixelDigitalCharge[Nhits_check][N_Thr_steps][MaxPixelInCluster];
  
  int Nbins_Mult = 40;
  double R_Mult[2];
  R_Mult[0] = -0.5;
  R_Mult[1] = Nbins_Mult + 0.5;
  Nbins_Mult++;
  TH1F* hmultiplicity[N_Thr_steps];
  
  int Nbins_SeedCharge = 100;
  double R_SeedCharge[2];
  R_SeedCharge[0] = 0.0;
  R_SeedCharge[1] = 5.0; //kelec
  TH1F* h_SeedCharge[N_Thr_steps];
  
  const int NMults(6);
  double Delta_factor = 1.5;
  double Fit_factor   = 3.0;
  int Nbins_DeltaPosition    = 100;
  double ThePitch = TMath::Max(GetPitchU(),GetPitchV());
  double R_DeltaPositionU[2];
  R_DeltaPositionU[0] = -Delta_factor*ThePitch;
  R_DeltaPositionU[1] = +Delta_factor*ThePitch;
  double R_DeltaPositionV[2];
  R_DeltaPositionV[0] = -Delta_factor*ThePitch;
  R_DeltaPositionV[1] = +Delta_factor*ThePitch;
  int Nbins_MyMult = NMults;
  double R_MyMult[2];
  R_MyMult[0] = 1.0 - 0.5;
  R_MyMult[1] = NMults + 0.5;
  TH1F* h_DeltaAnalogPositionU[N_Thr_steps];
  TH1F* h_DeltaAnalogPositionV[N_Thr_steps];
  TF1*  f_DeltaAnalogPositionU[N_Thr_steps];
  TF1*  f_DeltaAnalogPositionV[N_Thr_steps];
  TH2F* h_DeltaAnalogPositionUV[N_Thr_steps];
  TH1F* h_DeltaDigitalPositionU[N_Thr_steps];
  TH1F* h_DeltaDigitalPositionV[N_Thr_steps];
  TF1*  f_DeltaDigitalPositionU[N_Thr_steps];
  TF1*  f_DeltaDigitalPositionV[N_Thr_steps];
  TH2F* h_DeltaDigitalPositionUV[N_Thr_steps];
  TH1F* h_DeltaDigitalPositionU_vs_Mult[N_Thr_steps][NMults];
  TH1F* h_DeltaDigitalPositionV_vs_Mult[N_Thr_steps][NMults];
  TF1*  f_DeltaDigitalPositionU_vs_Mult[N_Thr_steps][NMults];
  TF1*  f_DeltaDigitalPositionV_vs_Mult[N_Thr_steps][NMults];
  TH1F* h_MeanDigitalDeltaU_vs_Mult[N_Thr_steps];
  TH1F* h_MeanDigitalDeltaV_vs_Mult[N_Thr_steps];
  TH1F* h_RMSDigitalDeltaU_vs_Mult[N_Thr_steps];
  TH1F* h_RMSDigitalDeltaV_vs_Mult[N_Thr_steps];
  
  int Nbins_TrkDistToDiode = 100;
  double R_TrkDistToDiode[2];
  R_TrkDistToDiode[0] = 0.0;
  R_TrkDistToDiode[1] = 1.5*TMath::Max(GetPitchU(),GetPitchV());
  TH1F* h_TrkDistToDiode_Detected[N_Thr_steps];
  TH1F* h_TrkDistToDiode_NonDetected[N_Thr_steps];
  TH1F* h_Effic_vs_TrkDistToDiode[N_Thr_steps];
  TH1F* h_TrkDistToDiode_Detected_vs_Mult[N_Thr_steps][NMults];
  
  int Nbins_TrkDistToDiode2D = 70;
  double R_TrkDistToDiode2DU[2];
  R_TrkDistToDiode2DU[0] = 0.0;
  R_TrkDistToDiode2DU[1] = 2.0*GetPitchU();
  double R_TrkDistToDiode2DV[2];
  R_TrkDistToDiode2DV[0] = 0.0;
  R_TrkDistToDiode2DV[1] = 2.0*GetPitchV();
  TH2F* h_TrkDistToDiode2D_Detected[N_Thr_steps];
  TH2F* h_TrkDistToDiode2D_NonDetected[N_Thr_steps];
  TH2F* h_TrkDistToDiode2D_Detected_vs_Mult[N_Thr_steps][NMults];
  
  TEllipse* Diodes[4];
  double rell_u = (R_TrkDistToDiode2DU[1] - R_TrkDistToDiode2DU[0])*0.03;
  double rell_v = (R_TrkDistToDiode2DV[1] - R_TrkDistToDiode2DV[0])*0.03;
  int col,lin;
  double u,v;
  col = 0; lin = 0;
  ComputePixelPositionUV_FromColRow(col,lin,u,v);
  u = (u + 0.5 * GetNpixelsU() * GetPitchU())/(2.0 * GetPitchU());
  u = (u - int(u))*2.0*GetPitchU();
  v = (v + 0.5 * GetNpixelsV() * GetPitchV())/(2.0 * GetPitchV());
  v = (v - int(v))*2.0*GetPitchV();
  Diodes[0] = new TEllipse(u,v,rell_u,rell_v);
  Diodes[0]->SetLineColor(2);
  Diodes[0]->SetLineWidth(2);
  Diodes[0]->SetFillStyle(3000);
  col = 0; lin = 1;
  ComputePixelPositionUV_FromColRow(col,lin,u,v);
  u = (u + 0.5 * GetNpixelsU() * GetPitchU())/(2.0 * GetPitchU());
  u = (u - int(u))*2.0*GetPitchU();
  v = (v + 0.5 * GetNpixelsV() * GetPitchV())/(2.0 * GetPitchV());
  v = (v - int(v))*2.0*GetPitchV();
  Diodes[1] = new TEllipse(u,v,rell_u,rell_v);
  Diodes[1]->SetLineColor(2);
  Diodes[1]->SetLineWidth(2);
  Diodes[1]->SetFillStyle(3000);
  col = 1; lin = 0;
  ComputePixelPositionUV_FromColRow(col,lin,u,v);
  u = (u + 0.5 * GetNpixelsU() * GetPitchU())/(2.0 * GetPitchU());
  u = (u - int(u))*2.0*GetPitchU();
  v = (v + 0.5 * GetNpixelsV() * GetPitchV())/(2.0 * GetPitchV());
  v = (v - int(v))*2.0*GetPitchV();
  Diodes[2] = new TEllipse(u,v,rell_u,rell_v);
  Diodes[2]->SetLineColor(2);
  Diodes[2]->SetLineWidth(2);
  Diodes[2]->SetFillStyle(3000);
  col = 1; lin = 1;
  ComputePixelPositionUV_FromColRow(col,lin,u,v);
  u = (u + 0.5 * GetNpixelsU() * GetPitchU())/(2.0 * GetPitchU());
  u = (u - int(u))*2.0*GetPitchU();
  v = (v + 0.5 * GetNpixelsV() * GetPitchV())/(2.0 * GetPitchV());
  v = (v - int(v))*2.0*GetPitchV();
  Diodes[3] = new TEllipse(u,v,rell_u,rell_v);
  Diodes[3]->SetLineColor(2);
  Diodes[3]->SetLineWidth(2);
  Diodes[3]->SetFillStyle(3000);
  
  
  TString HistName,HistTitle;
  
  HistName  = TString("h_effic_vs_Thr");
  HistTitle = TString("Detection efficiency vs threshold");
  TH1F* h_effic_vs_Thr = new TH1F(HistName.Data(),HistTitle.Data(),
				  N_Thr_steps,R_Thr[0],R_Thr[1]);
  h_effic_vs_Thr->SetXTitle("Threshold (#timesNoise)");
  h_effic_vs_Thr->GetXaxis()->CenterTitle(true);
  h_effic_vs_Thr->SetYTitle("Efficiency (%)");
  h_effic_vs_Thr->GetYaxis()->CenterTitle(true);
  h_effic_vs_Thr->SetLineColor(1);
  h_effic_vs_Thr->SetLineWidth(2);
  h_effic_vs_Thr->SetLineStyle(1);
  h_effic_vs_Thr->GetXaxis()->SetTitleSize(TheSize);
  h_effic_vs_Thr->GetXaxis()->SetLabelSize(TheSize);
  h_effic_vs_Thr->GetYaxis()->SetTitleSize(TheSize);
  h_effic_vs_Thr->GetYaxis()->SetLabelSize(TheSize);
  h_effic_vs_Thr->SetStats(false);
  
  HistName  = TString("h_AnalogResolutionU_vs_Thr");
  HistTitle = TString("Analog resolution vs threshold");
  TH1F* h_AnalogResolutionU_vs_Thr = new TH1F(HistName.Data(),HistTitle.Data(),
				              N_Thr_steps,R_Thr[0],R_Thr[1]);
  h_AnalogResolutionU_vs_Thr->SetXTitle("Threshold (#timesNoise)");
  h_AnalogResolutionU_vs_Thr->GetXaxis()->CenterTitle(true);
  h_AnalogResolutionU_vs_Thr->SetYTitle("Analog resolution (#mum)");
  h_AnalogResolutionU_vs_Thr->GetYaxis()->CenterTitle(true);
  h_AnalogResolutionU_vs_Thr->SetLineColor(kRed);
  h_AnalogResolutionU_vs_Thr->SetLineWidth(2);
  h_AnalogResolutionU_vs_Thr->SetLineStyle(1);
  h_AnalogResolutionU_vs_Thr->GetXaxis()->SetTitleSize(TheSize);
  h_AnalogResolutionU_vs_Thr->GetXaxis()->SetLabelSize(TheSize);
  h_AnalogResolutionU_vs_Thr->GetYaxis()->SetTitleSize(TheSize);
  h_AnalogResolutionU_vs_Thr->GetYaxis()->SetLabelSize(TheSize);
  h_AnalogResolutionU_vs_Thr->SetStats(false);
  
  HistName  = TString("h_AnalogResolutionV_vs_Thr");
  HistTitle = TString("Analog resolution vs threshold");
  TH1F* h_AnalogResolutionV_vs_Thr = new TH1F(HistName.Data(),HistTitle.Data(),
				              N_Thr_steps,R_Thr[0],R_Thr[1]);
  h_AnalogResolutionV_vs_Thr->SetXTitle("Threshold (#timesNoise)");
  h_AnalogResolutionV_vs_Thr->GetXaxis()->CenterTitle(true);
  h_AnalogResolutionV_vs_Thr->SetYTitle("Analog resolution (#mum)");
  h_AnalogResolutionV_vs_Thr->GetYaxis()->CenterTitle(true);
  h_AnalogResolutionV_vs_Thr->SetLineColor(kRed);
  h_AnalogResolutionV_vs_Thr->SetLineWidth(2);
  h_AnalogResolutionV_vs_Thr->SetLineStyle(2);
  h_AnalogResolutionV_vs_Thr->GetXaxis()->SetTitleSize(TheSize);
  h_AnalogResolutionV_vs_Thr->GetXaxis()->SetLabelSize(TheSize);
  h_AnalogResolutionV_vs_Thr->GetYaxis()->SetTitleSize(TheSize);
  h_AnalogResolutionV_vs_Thr->GetYaxis()->SetLabelSize(TheSize);
  h_AnalogResolutionV_vs_Thr->SetStats(false);
  
  HistName  = TString("h_DigitalResolutionU_vs_Thr");
  HistTitle = TString("Digital resolution vs threshold");
  TH1F* h_DigitalResolutionU_vs_Thr = new TH1F(HistName.Data(),HistTitle.Data(),
				               N_Thr_steps,R_Thr[0],R_Thr[1]);
  h_DigitalResolutionU_vs_Thr->SetXTitle("Threshold (#timesNoise)");
  h_DigitalResolutionU_vs_Thr->GetXaxis()->CenterTitle(true);
  h_DigitalResolutionU_vs_Thr->SetYTitle("Ditital resolution (#mum)");
  h_DigitalResolutionU_vs_Thr->GetYaxis()->CenterTitle(true);
  h_DigitalResolutionU_vs_Thr->SetLineColor(kRed);
  h_DigitalResolutionU_vs_Thr->SetLineWidth(2);
  h_DigitalResolutionU_vs_Thr->SetLineStyle(1);
  h_DigitalResolutionU_vs_Thr->GetXaxis()->SetTitleSize(TheSize);
  h_DigitalResolutionU_vs_Thr->GetXaxis()->SetLabelSize(TheSize);
  h_DigitalResolutionU_vs_Thr->GetYaxis()->SetTitleSize(TheSize);
  h_DigitalResolutionU_vs_Thr->GetYaxis()->SetLabelSize(TheSize);
  h_DigitalResolutionU_vs_Thr->SetStats(false);
  
  HistName  = TString("h_DigitalResolutionV_vs_Thr");
  HistTitle = TString("Digital resolution vs threshold");
  TH1F* h_DigitalResolutionV_vs_Thr = new TH1F(HistName.Data(),HistTitle.Data(),
				               N_Thr_steps,R_Thr[0],R_Thr[1]);
  h_DigitalResolutionV_vs_Thr->SetXTitle("Threshold (#timesNoise)");
  h_DigitalResolutionV_vs_Thr->GetXaxis()->CenterTitle(true);
  h_DigitalResolutionV_vs_Thr->SetYTitle("Digital resolution V (#mum)");
  h_DigitalResolutionV_vs_Thr->GetYaxis()->CenterTitle(true);
  h_DigitalResolutionV_vs_Thr->SetLineColor(kRed);
  h_DigitalResolutionV_vs_Thr->SetLineWidth(2);
  h_DigitalResolutionV_vs_Thr->SetLineStyle(2);
  h_DigitalResolutionV_vs_Thr->GetXaxis()->SetTitleSize(TheSize);
  h_DigitalResolutionV_vs_Thr->GetXaxis()->SetLabelSize(TheSize);
  h_DigitalResolutionV_vs_Thr->GetYaxis()->SetTitleSize(TheSize);
  h_DigitalResolutionV_vs_Thr->GetYaxis()->SetLabelSize(TheSize);
  h_DigitalResolutionV_vs_Thr->SetStats(false);

  char ytitle[100];
  char ytitle2[100];
  for(int ithr=0;ithr<N_Thr_steps;ithr++) {
    double Threshold_Times_Noise = R_Thr[0] + (ithr + 0.5)*(R_Thr[1] - R_Thr[0])/N_Thr_steps;
    sprintf(ytitle,"%.1f",Threshold_Times_Noise);
    
    HistName  = TString("hmultiplicity_Thr") + long(ithr+1);
    HistTitle = TString("cluster multiplicity for Thr = ") + TString(ytitle) + TString(" #times Noise");
    hmultiplicity[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                   Nbins_Mult,R_Mult[0],R_Mult[1]);
    hmultiplicity[ithr]->SetXTitle("# pixels / cluster");
    hmultiplicity[ithr]->GetXaxis()->CenterTitle(true);
    hmultiplicity[ithr]->SetYTitle("Hits");
    hmultiplicity[ithr]->GetYaxis()->CenterTitle(true);
    hmultiplicity[ithr]->SetLineColor(4);
    hmultiplicity[ithr]->SetLineWidth(2);
    hmultiplicity[ithr]->SetMinimum(1.0e-3);
    hmultiplicity[ithr]->GetXaxis()->SetTitleSize(TheSize);
    hmultiplicity[ithr]->GetXaxis()->SetLabelSize(TheSize);
    hmultiplicity[ithr]->GetYaxis()->SetTitleSize(TheSize);
    hmultiplicity[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_SeedCharge_Thr") + long(ithr+1);
    HistTitle = TString("Seed pixel charge for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_SeedCharge[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                  Nbins_SeedCharge,R_SeedCharge[0],R_SeedCharge[1]);
    h_SeedCharge[ithr]->SetXTitle("Seed pixel charge (#times10^{3} electrons)");
    h_SeedCharge[ithr]->GetXaxis()->CenterTitle(true);
    h_SeedCharge[ithr]->SetYTitle("Hits");
    h_SeedCharge[ithr]->GetYaxis()->CenterTitle(true);
    h_SeedCharge[ithr]->SetLineColor(4);
    h_SeedCharge[ithr]->SetLineWidth(2);
    h_SeedCharge[ithr]->SetMinimum(1.0e-3);
    h_SeedCharge[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_SeedCharge[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_SeedCharge[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_SeedCharge[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_DeltaAnalogPositionU_Thr") + long(ithr+1);
    HistTitle = TString("Analog #DeltaU for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_DeltaAnalogPositionU[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                            Nbins_DeltaPosition,R_DeltaPositionU[0],R_DeltaPositionU[1]);
    h_DeltaAnalogPositionU[ithr]->SetXTitle("Analog #DeltaU (#mum)");
    h_DeltaAnalogPositionU[ithr]->GetXaxis()->CenterTitle(true);
    h_DeltaAnalogPositionU[ithr]->SetYTitle("Hits");
    h_DeltaAnalogPositionU[ithr]->GetYaxis()->CenterTitle(true);
    h_DeltaAnalogPositionU[ithr]->SetLineColor(kBlue);
    h_DeltaAnalogPositionU[ithr]->SetLineWidth(2);
    h_DeltaAnalogPositionU[ithr]->SetMinimum(1.0e-3);
    h_DeltaAnalogPositionU[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionU[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_DeltaAnalogPositionU[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionU[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_DeltaAnalogPositionV_Thr") + long(ithr+1);
    HistTitle = TString("Analog #DeltaV for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_DeltaAnalogPositionV[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                            Nbins_DeltaPosition,R_DeltaPositionV[0],R_DeltaPositionV[1]);
    h_DeltaAnalogPositionV[ithr]->SetXTitle("Analog #DeltaV (#mum)");
    h_DeltaAnalogPositionV[ithr]->GetXaxis()->CenterTitle(true);
    h_DeltaAnalogPositionV[ithr]->SetYTitle("Hits");
    h_DeltaAnalogPositionV[ithr]->GetYaxis()->CenterTitle(true);
    h_DeltaAnalogPositionV[ithr]->SetLineColor(kBlue);
    h_DeltaAnalogPositionV[ithr]->SetLineWidth(2);
    h_DeltaAnalogPositionV[ithr]->SetMinimum(1.0e-3);
    h_DeltaAnalogPositionV[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionV[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_DeltaAnalogPositionV[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionV[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_DeltaAnalogPositionUV_Thr") + long(ithr+1);
    HistTitle = TString("Analog #DeltaU vs #DeltaV for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_DeltaAnalogPositionUV[ithr] = new TH2F(HistName.Data(),HistTitle.Data(),
					     Nbins_DeltaPosition,R_DeltaPositionU[0],R_DeltaPositionU[1],
                                             Nbins_DeltaPosition,R_DeltaPositionV[0],R_DeltaPositionV[1]);
    h_DeltaAnalogPositionUV[ithr]->SetXTitle("Analog #DeltaU (#mum)");
    h_DeltaAnalogPositionUV[ithr]->GetXaxis()->CenterTitle(true);
    h_DeltaAnalogPositionUV[ithr]->SetYTitle("Analog #DeltaV (#mum)");
    h_DeltaAnalogPositionUV[ithr]->GetYaxis()->CenterTitle(true);
    h_DeltaAnalogPositionUV[ithr]->SetZTitle("Hits");
    h_DeltaAnalogPositionUV[ithr]->GetZaxis()->CenterTitle(true);
    h_DeltaAnalogPositionUV[ithr]->SetLineColor(kBlue);
    h_DeltaAnalogPositionUV[ithr]->SetLineWidth(2);
    h_DeltaAnalogPositionUV[ithr]->SetMinimum(1.0e-3);
    h_DeltaAnalogPositionUV[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionUV[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_DeltaAnalogPositionUV[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionUV[ithr]->GetYaxis()->SetLabelSize(TheSize);
    h_DeltaAnalogPositionUV[ithr]->GetZaxis()->SetTitleSize(TheSize);
    h_DeltaAnalogPositionUV[ithr]->GetZaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_DeltaDigitalPositionU_Thr") + long(ithr+1);
    HistTitle = TString("Digital #DeltaU for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_DeltaDigitalPositionU[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                             Nbins_DeltaPosition,R_DeltaPositionU[0],R_DeltaPositionU[1]);
    h_DeltaDigitalPositionU[ithr]->SetXTitle("Digital #DeltaU (#mum)");
    h_DeltaDigitalPositionU[ithr]->GetXaxis()->CenterTitle(true);
    h_DeltaDigitalPositionU[ithr]->SetYTitle("Hits");
    h_DeltaDigitalPositionU[ithr]->GetYaxis()->CenterTitle(true);
    h_DeltaDigitalPositionU[ithr]->SetLineColor(kBlue);
    h_DeltaDigitalPositionU[ithr]->SetLineWidth(2);
    h_DeltaDigitalPositionU[ithr]->SetMinimum(1.0e-3);
    h_DeltaDigitalPositionU[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionU[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_DeltaDigitalPositionU[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionU[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_DeltaDigitalPositionV_Thr") + long(ithr+1);
    HistTitle = TString("Digital #DeltaV for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_DeltaDigitalPositionV[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                             Nbins_DeltaPosition,R_DeltaPositionV[0],R_DeltaPositionV[1]);
    h_DeltaDigitalPositionV[ithr]->SetXTitle("Digital #DeltaV (#mum)");
    h_DeltaDigitalPositionV[ithr]->GetXaxis()->CenterTitle(true);
    h_DeltaDigitalPositionV[ithr]->SetYTitle("Hits");
    h_DeltaDigitalPositionV[ithr]->GetYaxis()->CenterTitle(true);
    h_DeltaDigitalPositionV[ithr]->SetLineColor(kBlue);
    h_DeltaDigitalPositionV[ithr]->SetLineWidth(2);
    h_DeltaDigitalPositionV[ithr]->SetMinimum(1.0e-3);
    h_DeltaDigitalPositionV[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionV[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_DeltaDigitalPositionV[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionV[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_DeltaDigitalPositionUV_Thr") + long(ithr+1);
    HistTitle = TString("Digital #DeltaU vs #DeltaV for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_DeltaDigitalPositionUV[ithr] = new TH2F(HistName.Data(),HistTitle.Data(),
		 			      Nbins_DeltaPosition,R_DeltaPositionU[0],R_DeltaPositionU[1],
                                              Nbins_DeltaPosition,R_DeltaPositionV[0],R_DeltaPositionV[1]);
    h_DeltaDigitalPositionUV[ithr]->SetXTitle("Digital #DeltaU (#mum)");
    h_DeltaDigitalPositionUV[ithr]->GetXaxis()->CenterTitle(true);
    h_DeltaDigitalPositionUV[ithr]->SetYTitle("Digital #DeltaV (#mum)");
    h_DeltaDigitalPositionUV[ithr]->GetYaxis()->CenterTitle(true);
    h_DeltaDigitalPositionUV[ithr]->SetZTitle("Hits");
    h_DeltaDigitalPositionUV[ithr]->GetZaxis()->CenterTitle(true);
    h_DeltaDigitalPositionUV[ithr]->SetLineColor(kBlue);
    h_DeltaDigitalPositionUV[ithr]->SetLineWidth(2);
    h_DeltaDigitalPositionUV[ithr]->SetMinimum(1.0e-3);
    h_DeltaDigitalPositionUV[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionUV[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_DeltaDigitalPositionUV[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionUV[ithr]->GetYaxis()->SetLabelSize(TheSize);
    h_DeltaDigitalPositionUV[ithr]->GetZaxis()->SetTitleSize(TheSize);
    h_DeltaDigitalPositionUV[ithr]->GetZaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_MeanDigitalDeltaU_vs_Mult_Thr") + long(ithr+1);
    HistTitle = TString("mean #DeltaU vs mult for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_MeanDigitalDeltaU_vs_Mult[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                Nbins_MyMult,R_MyMult[0],R_MyMult[1]);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->SetXTitle("# pixels / cluster");
    h_MeanDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->CenterTitle(true);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->SetYTitle("#mu(#DeltaU/V) (#mum)");
    h_MeanDigitalDeltaU_vs_Mult[ithr]->GetYaxis()->CenterTitle(true);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->SetLineColor(kBlue);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->SetLineWidth(2);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_MeanDigitalDeltaU_vs_Mult[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_MeanDigitalDeltaV_vs_Mult_Thr") + long(ithr+1);
    HistTitle = TString("mean #DeltaV vs mult for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_MeanDigitalDeltaV_vs_Mult[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                Nbins_MyMult,R_MyMult[0],R_MyMult[1]);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->SetXTitle("# pixels / cluster");
    h_MeanDigitalDeltaV_vs_Mult[ithr]->GetXaxis()->CenterTitle(true);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->SetYTitle("#mu(#DeltaU/V) (#mum)");
    h_MeanDigitalDeltaV_vs_Mult[ithr]->GetYaxis()->CenterTitle(true);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->SetLineColor(kRed);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->SetLineWidth(2);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_MeanDigitalDeltaV_vs_Mult[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_RMSDigitalDeltaU_vs_Mult_Thr") + long(ithr+1);
    HistTitle = TString("rms #DeltaU vs mult for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_RMSDigitalDeltaU_vs_Mult[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                Nbins_MyMult,R_MyMult[0],R_MyMult[1]);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->SetXTitle("# pixels / cluster");
    h_RMSDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->CenterTitle(true);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->SetYTitle("#sigma(#DeltaU/V) (#mum)");
    h_RMSDigitalDeltaU_vs_Mult[ithr]->GetYaxis()->CenterTitle(true);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->SetLineColor(kBlue);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->SetLineWidth(2);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_RMSDigitalDeltaU_vs_Mult[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_RMSDigitalDeltaV_vs_Mult_Thr") + long(ithr+1);
    HistTitle = TString("rms #DeltaV vs mult for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_RMSDigitalDeltaV_vs_Mult[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                Nbins_MyMult,R_MyMult[0],R_MyMult[1]);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->SetXTitle("# pixels / cluster");
    h_RMSDigitalDeltaV_vs_Mult[ithr]->GetXaxis()->CenterTitle(true);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->SetYTitle("#sigma(#DeltaU/V) (#mum)");
    h_RMSDigitalDeltaV_vs_Mult[ithr]->GetYaxis()->CenterTitle(true);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->SetLineColor(kRed);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->SetLineWidth(2);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_RMSDigitalDeltaV_vs_Mult[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
  
    HistName  = TString("h_TrkDistToDiode_Detected_Thr") + long(ithr+1);
    HistTitle = TString("Trk distance to closest diode for all detected hits for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_TrkDistToDiode_Detected[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                               Nbins_TrkDistToDiode,R_TrkDistToDiode[0],R_TrkDistToDiode[1]);
    h_TrkDistToDiode_Detected[ithr]->SetXTitle("Trk dist to closest diode (#mum)");
    h_TrkDistToDiode_Detected[ithr]->GetXaxis()->CenterTitle(true);
    h_TrkDistToDiode_Detected[ithr]->SetYTitle("Hits");
    h_TrkDistToDiode_Detected[ithr]->GetYaxis()->CenterTitle(true);
    h_TrkDistToDiode_Detected[ithr]->SetLineColor(kBlue);
    h_TrkDistToDiode_Detected[ithr]->SetLineWidth(2);
    h_TrkDistToDiode_Detected[ithr]->SetMinimum(1.0e-3);
    h_TrkDistToDiode_Detected[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode_Detected[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode_Detected[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode_Detected[ithr]->GetYaxis()->SetLabelSize(TheSize);
    
    HistName  = TString("h_TrkDistToDiode_NonDetected_Thr") + long(ithr+1);
    HistTitle = TString("Trk distance to closest diode for non-detected hits for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_TrkDistToDiode_NonDetected[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                  Nbins_TrkDistToDiode,R_TrkDistToDiode[0],R_TrkDistToDiode[1]);
    h_TrkDistToDiode_NonDetected[ithr]->SetXTitle("Trk dist to closest diode (#mum)");
    h_TrkDistToDiode_NonDetected[ithr]->GetXaxis()->CenterTitle(true);
    h_TrkDistToDiode_NonDetected[ithr]->SetYTitle("Hits");
    h_TrkDistToDiode_NonDetected[ithr]->GetYaxis()->CenterTitle(true);
    h_TrkDistToDiode_NonDetected[ithr]->SetLineColor(kBlue);
    h_TrkDistToDiode_NonDetected[ithr]->SetLineWidth(2);
    h_TrkDistToDiode_NonDetected[ithr]->SetMinimum(1.0e-3);
    h_TrkDistToDiode_NonDetected[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode_NonDetected[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode_NonDetected[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode_NonDetected[ithr]->GetYaxis()->SetLabelSize(TheSize);
  
    HistName  = TString("h_Effic_vs_TrkDistToDiode_Thr") + long(ithr+1);
    HistTitle = TString("Efficiency vs Trk distance to closest diode for all detected hits for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_Effic_vs_TrkDistToDiode[ithr] = new TH1F(HistName.Data(),HistTitle.Data(),
                                               Nbins_TrkDistToDiode,R_TrkDistToDiode[0],R_TrkDistToDiode[1]);
    h_Effic_vs_TrkDistToDiode[ithr]->SetXTitle("Trk dist to closest diode (#mum)");
    h_Effic_vs_TrkDistToDiode[ithr]->GetXaxis()->CenterTitle(true);
    h_Effic_vs_TrkDistToDiode[ithr]->SetYTitle("Efficiency (%)");
    h_Effic_vs_TrkDistToDiode[ithr]->GetYaxis()->CenterTitle(true);
    h_Effic_vs_TrkDistToDiode[ithr]->SetLineColor(kBlue);
    h_Effic_vs_TrkDistToDiode[ithr]->SetLineWidth(2);
    h_Effic_vs_TrkDistToDiode[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_Effic_vs_TrkDistToDiode[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_Effic_vs_TrkDistToDiode[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_Effic_vs_TrkDistToDiode[ithr]->GetYaxis()->SetLabelSize(TheSize);
    h_Effic_vs_TrkDistToDiode[ithr]->SetStats(false);
  
    HistName  = TString("h_TrkDistToDiode2D_Detected_Thr") + long(ithr+1);
    HistTitle = TString("Trk U vs V position for all detected hits for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_TrkDistToDiode2D_Detected[ithr] = new TH2F(HistName.Data(),HistTitle.Data(),
		 			         Nbins_TrkDistToDiode2D,R_TrkDistToDiode2DU[0],R_TrkDistToDiode2DU[1],
                                                 Nbins_TrkDistToDiode2D,R_TrkDistToDiode2DV[0],R_TrkDistToDiode2DV[1]);
    h_TrkDistToDiode2D_Detected[ithr]->SetXTitle("Trk U position (#mum)");
    h_TrkDistToDiode2D_Detected[ithr]->GetXaxis()->CenterTitle(true);
    h_TrkDistToDiode2D_Detected[ithr]->SetYTitle("Trk V position (#mum)");
    h_TrkDistToDiode2D_Detected[ithr]->GetYaxis()->CenterTitle(true);
    h_TrkDistToDiode2D_Detected[ithr]->SetZTitle("Hits");
    h_TrkDistToDiode2D_Detected[ithr]->GetZaxis()->CenterTitle(true);
    h_TrkDistToDiode2D_Detected[ithr]->SetLineColor(kBlue);
    h_TrkDistToDiode2D_Detected[ithr]->SetLineWidth(2);
    h_TrkDistToDiode2D_Detected[ithr]->SetMinimum(1.0e-3);
    h_TrkDistToDiode2D_Detected[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode2D_Detected[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode2D_Detected[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode2D_Detected[ithr]->GetYaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode2D_Detected[ithr]->GetZaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode2D_Detected[ithr]->GetZaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode2D_Detected[ithr]->SetStats(false);
    
    HistName  = TString("h_TrkDistToDiode2D_NonDetected_Thr") + long(ithr+1);
    HistTitle = TString("Trk U vs V position for all detected hits for Thr = ") + TString(ytitle) + TString(" #times Noise");
    h_TrkDistToDiode2D_NonDetected[ithr] = new TH2F(HistName.Data(),HistTitle.Data(),
		 			            Nbins_TrkDistToDiode2D,R_TrkDistToDiode2DU[0],R_TrkDistToDiode2DU[1],
                                                    Nbins_TrkDistToDiode2D,R_TrkDistToDiode2DV[0],R_TrkDistToDiode2DV[1]);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetXTitle("Trk U position (#mum)");
    h_TrkDistToDiode2D_NonDetected[ithr]->GetXaxis()->CenterTitle(true);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetYTitle("Trk V position (#mum)");
    h_TrkDistToDiode2D_NonDetected[ithr]->GetYaxis()->CenterTitle(true);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetZTitle("Hits");
    h_TrkDistToDiode2D_NonDetected[ithr]->GetZaxis()->CenterTitle(true);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetLineColor(kBlue);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetLineWidth(2);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetMinimum(1.0e-3);
    h_TrkDistToDiode2D_NonDetected[ithr]->GetXaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode2D_NonDetected[ithr]->GetXaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode2D_NonDetected[ithr]->GetYaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode2D_NonDetected[ithr]->GetYaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode2D_NonDetected[ithr]->GetZaxis()->SetTitleSize(TheSize);
    h_TrkDistToDiode2D_NonDetected[ithr]->GetZaxis()->SetLabelSize(TheSize);
    h_TrkDistToDiode2D_NonDetected[ithr]->SetStats(false);
    
    for(int imult=0;imult<NMults;imult++) {
      HistName  = TString("h_DeltaDigitalPositionU_Thr") + long(ithr+1) + TString("_mult") + long(imult+1);
      if(imult+1 < NMults) HistTitle = TString("Digital #DeltaU mult =  ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      else                 HistTitle = TString("Digital #DeltaU mult >= ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                              Nbins_DeltaPosition,R_DeltaPositionU[0],R_DeltaPositionU[1]);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetXTitle("Digital #DeltaU (#mum)");
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetXaxis()->CenterTitle(true);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetYTitle("Hits");
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetYaxis()->CenterTitle(true);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetLineColor(kBlue);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetLineWidth(2);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetMinimum(1.0e-3);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetXaxis()->SetTitleSize(TheSize);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetXaxis()->SetLabelSize(TheSize);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetYaxis()->SetTitleSize(TheSize);
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetYaxis()->SetLabelSize(TheSize);
      
      HistName  = TString("h_DeltaDigitalPositionV_Thr") + long(ithr+1) + TString("_mult") + long(imult+1);
      if(imult+1 < NMults) HistTitle = TString("Digital #DeltaV mult =  ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      else                 HistTitle = TString("Digital #DeltaV mult >= ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                              Nbins_DeltaPosition,R_DeltaPositionU[0],R_DeltaPositionU[1]);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetXTitle("Digital #DeltaV (#mum)");
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetXaxis()->CenterTitle(true);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetYTitle("Hits");
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetYaxis()->CenterTitle(true);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetLineColor(kBlue);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetLineWidth(2);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetMinimum(1.0e-3);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetXaxis()->SetTitleSize(TheSize);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetXaxis()->SetLabelSize(TheSize);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetYaxis()->SetTitleSize(TheSize);
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetYaxis()->SetLabelSize(TheSize);
      
      HistName  = TString("h_TrkDistToDiode_Detected_Thr") + long(ithr+1) + TString("_mult") + long(imult+1);
      if(imult+1 < NMults) HistTitle = TString("Trk distance to closest diode for detected hits for mult =  ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      else                 HistTitle = TString("Trk distance to closest diode for detected hits for mult >= ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult] = new TH1F(HistName.Data(),HistTitle.Data(),
                                                                Nbins_TrkDistToDiode,R_TrkDistToDiode[0],R_TrkDistToDiode[1]);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->SetXTitle("Trk dist to closest diode (#mum)");
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->GetXaxis()->CenterTitle(true);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->SetYTitle("Hits");
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->GetYaxis()->CenterTitle(true);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->SetLineColor(kBlue);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->SetLineWidth(2);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->SetMinimum(1.0e-3);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->GetXaxis()->SetTitleSize(TheSize);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->GetXaxis()->SetLabelSize(TheSize);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->GetYaxis()->SetTitleSize(TheSize);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->GetYaxis()->SetLabelSize(TheSize);
      
      HistName  = TString("h_TrkDistToDiode2D_Detected_Thr") + long(ithr+1) + TString("_mult") + long(imult+1);
      if(imult+1 < NMults) HistTitle = TString("Trk U vs V position for detected hits for mult =  ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      else                 HistTitle = TString("Trk U vs V position for detected hits for mult >= ") + long(imult+1) + (" and for Thr = ") + TString(ytitle) + TString(" #times Noise");
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult] = new TH2F(HistName.Data(),HistTitle.Data(),
		 			                          Nbins_TrkDistToDiode2D,R_TrkDistToDiode2DU[0],R_TrkDistToDiode2DU[1],
                                                                  Nbins_TrkDistToDiode2D,R_TrkDistToDiode2DV[0],R_TrkDistToDiode2DV[1]);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetXTitle("Trk U position (#mum)");
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetXaxis()->CenterTitle(true);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetYTitle("Trk V position (#mum)");
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetYaxis()->CenterTitle(true);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetZTitle("Hits");
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetZaxis()->CenterTitle(true);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetLineColor(kBlue);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetLineWidth(2);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetMinimum(1.0e-3);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetXaxis()->SetTitleSize(TheSize);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetXaxis()->SetLabelSize(TheSize);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetYaxis()->SetTitleSize(TheSize);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetYaxis()->SetLabelSize(TheSize);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetZaxis()->SetTitleSize(TheSize);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->GetZaxis()->SetLabelSize(TheSize);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->SetStats(false);
    }
    
    Double_t Discri_Threshold_electrons = Threshold_Times_Noise*NoiseElectrons;
    sprintf(ytitle2,"%.3f",Discri_Threshold_electrons*1.0e-3);
    int Ndet_evts     = 0;
    NhitsToDraw[ithr] = 0;
    
    cout << endl;
    cout << "Doing simulation for Threshold = " << Threshold_Times_Noise << " x Noise (" << ithr+1 << " out of " << N_Thr_steps << ")" << endl;
    for(Int_t Nev = 0; Nev <Nevents; Nev++){
      if(!((Nev+1)%fPrintFreq)) {
        cout << Nev+1 << " simulated events!!!  ";
        fWatch.Print();
        fWatch.Continue();
      }
      
      Double_t RandU     = r3->Rndm()*(RU[1] - RU[0]) + RU[0];
      Double_t RandV     = r3->Rndm()*(RV[1] - RV[0]) + RV[0];

      Double_t RandTheta = Trk_thetaLoc*TMath::Pi()/180.0;
      Double_t RandPhi   = r3->Rndm()*2.0*TMath::Pi();

      Double_t InputUpos    = RandU;
      Double_t InputVpos    = RandV;
      Double_t totalUlength = EffectiveEpi*TMath::Tan(RandTheta)*TMath::Cos(RandPhi);
      Double_t totalVlength = EffectiveEpi*TMath::Tan(RandTheta)*TMath::Sin(RandPhi);
      Double_t OutputUpos   = InputUpos + totalUlength;
      Double_t OutputVpos   = InputVpos + totalVlength;
      
      Double_t AveUpos      = 0.5*(InputUpos + OutputUpos);
      Double_t AveVpos      = 0.5*(InputVpos + OutputVpos);

      double trk_dist_to_diode = GetTrkDistToDiode(AveUpos,AveVpos);
      TVector2 ShiftedPosition = GetTrkShiftedPosition(AveUpos,AveVpos);
      
      //pixel with deposited charge list
      fPixelMap.clear();
      fAnalogChargeMap.clear();
      fDigitalChargeMap.clear();
      Transport_Noise_Digitisation(SegmentSize,
				   InputUpos,InputVpos,
				   OutputUpos,OutputVpos,
				   -1.0,
				   RandTheta,
				   RandPhi,
				   fPixelMap,fAnalogChargeMap,fDigitalChargeMap,
				   Discri_Threshold_electrons);

      if(verbose) {
        cout << endl;
        cout << "Event " << Nev+1 << " begin:"<< endl;
        cout << " - (theta,phi)  = (" << RandTheta*180.0/TMath::Pi() << "," << RandPhi*180.0/TMath::Pi() << ") deg" << endl;
        cout << " - In-Pos(u,v)  = (" << InputUpos  << "," << InputVpos  << ") um" << endl;
        cout << " - Out-Pos(u,v) = (" << OutputUpos << "," << OutputVpos << ") um" << endl;
      }

      
      if(verbose) {
        cout << "Printing out the " << fPixelMap.size() << " found pixels:" << endl;
      }
      
      double R_ClusterCheck_Col[2];
      double R_ClusterCheck_Row[2];
      R_ClusterCheck_Col[0] = +1.0e+20;
      R_ClusterCheck_Col[1] = -1.0e+20;
      R_ClusterCheck_Row[0] = +1.0e+20;
      R_ClusterCheck_Row[1] = -1.0e+20;
      
      //NhitsToDraw[ithr] = 0;
      //TH2F* h_ClusterShape[Nhits_check][N_Thr_steps];
      
      Int_t  totalMult            = 0;
      double Seed_Charge          = -1.0e+20;
      double Position_analog_U    = 0.0;
      double Position_analog_V    = 0.0;
      double SumOfWeights_analog  = 0.0;
      double Position_digital_U   = 0.0;
      double Position_digital_V   = 0.0;
      bool   IsDetected_analog    = false;
      for(Int_t ipix = 0; ipix < int(fPixelMap.size()); ipix++){
        Int_t col,row;
        Double_t u,v;
	GetXYPixelNumber(col,row,fPixelMap[ipix]);
        ComputePixelPositionUV_FromColRow(col,row,u,v);
      
	if(R_ClusterCheck_Col[0] > col) R_ClusterCheck_Col[0] = col;
	if(R_ClusterCheck_Col[1] < col) R_ClusterCheck_Col[1] = col;
	if(R_ClusterCheck_Row[0] > row) R_ClusterCheck_Row[0] = row;
	if(R_ClusterCheck_Row[1] < row) R_ClusterCheck_Row[1] = row;
	
        if(verbose) {
	  cout << "pixel " << ipix << ", " 
	       << "(col,row) = (" << col << "," << row << "), " 
	       << "(u,v) = (" << u << "," << v << ") um, "
	       << "Analog Charge = " << fAnalogChargeMap[ipix] << " elec, "
	       << "Digital Charge = " << fDigitalChargeMap[ipix]
	       << endl;
        }
        
        //Cluster position analog
        Position_analog_U   += u*fAnalogChargeMap[ipix];
	Position_analog_V   += v*fAnalogChargeMap[ipix];
	SumOfWeights_analog +=   fAnalogChargeMap[ipix];
	
	//Cluster position analog
	if(fDigitalChargeMap[ipix] > 0) {
          Position_digital_U += u;
	  Position_digital_V += v;
	}
        
        totalMult += fDigitalChargeMap[ipix];
	if(Seed_Charge < fAnalogChargeMap[ipix]) Seed_Charge = fAnalogChargeMap[ipix];
      }
      
      if(NhitsToDraw[ithr] < Nhits_check && NhitsToDraw[ithr] < Nevents) {
	R_ClusterCheck_Col[0] -= 1;
	R_ClusterCheck_Col[1] += 1;
	R_ClusterCheck_Row[0] -= 1;
	R_ClusterCheck_Row[1] += 1;
	R_ClusterCheck_Col[0] -= 0.5;
	R_ClusterCheck_Col[1] += 0.5;
	R_ClusterCheck_Row[0] -= 0.5;
	R_ClusterCheck_Row[1] += 0.5;
	
	int Nbins_ClusterCheck_Col = int(R_ClusterCheck_Col[1] - R_ClusterCheck_Col[0]);
	int Nbins_ClusterCheck_Row = int(R_ClusterCheck_Row[1] - R_ClusterCheck_Row[0]);
	
	
	HistName  = TString("h_ClusterShape_Analog_hit") + long(NhitsToDraw[ithr]+1) + TString("_Thr") + long(ithr+1);
	HistTitle = TString("Analog 2D cluster for hit ") + long(NhitsToDraw[ithr]+1) + (" for Thr = ") + TString(ytitle) + TString(" #times Noise (") + TString(ytitle2) + TString("#times10^{3} elec)");
	h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr] = new TH2F(HistName.Data(),HistTitle.Data(),
								  Nbins_ClusterCheck_Col,R_ClusterCheck_Col[0],R_ClusterCheck_Col[1],
							          Nbins_ClusterCheck_Row,R_ClusterCheck_Row[0],R_ClusterCheck_Row[1]);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetXTitle("col");
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetXaxis()->CenterTitle(true);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetYTitle("row");
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetYaxis()->CenterTitle(true);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetZTitle("Analog Charge (#times10^{3} electrons)");
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetZaxis()->CenterTitle(true);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetLineColor(kBlue);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetLineWidth(2);
        //h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetMinimum(1.0e-3);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetXaxis()->SetTitleSize(TheSize);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetXaxis()->SetLabelSize(TheSize);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetYaxis()->SetTitleSize(TheSize);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetYaxis()->SetLabelSize(TheSize);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetZaxis()->SetTitleSize(TheSize);
        h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetZaxis()->SetLabelSize(TheSize);
	h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetStats(false);
	
	HistName  = TString("h_ClusterShape_Digital_hit") + long(NhitsToDraw[ithr]+1) + TString("_Thr") + long(ithr+1);
        HistTitle = TString("Digital 2D cluster for hit ") + long(NhitsToDraw[ithr]+1) + (" for Thr = ") + TString(ytitle) + TString(" #times Noise (") + TString(ytitle2) + TString("#times10^{3} elec)");
	h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr] = new TH2F(HistName.Data(),HistTitle.Data(),
								   Nbins_ClusterCheck_Col,R_ClusterCheck_Col[0],R_ClusterCheck_Col[1],
							           Nbins_ClusterCheck_Row,R_ClusterCheck_Row[0],R_ClusterCheck_Row[1]);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetXTitle("col");
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetXaxis()->CenterTitle(true);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetYTitle("row");
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetYaxis()->CenterTitle(true);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetZTitle("");
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetZaxis()->CenterTitle(true);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetLineColor(kBlue);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetLineWidth(2);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetMinimum(1.0e-3);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetXaxis()->SetTitleSize(TheSize);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetXaxis()->SetLabelSize(TheSize);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetYaxis()->SetTitleSize(TheSize);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetYaxis()->SetLabelSize(TheSize);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetZaxis()->SetTitleSize(TheSize);
        h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->GetZaxis()->SetLabelSize(TheSize);
	h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetStats(false);
	
	double rell_ppp   = 0.015*TMath::Max(Nbins_ClusterCheck_Col,Nbins_ClusterCheck_Row);
	double rell_pppUV = 0.015*TMath::Max(Nbins_ClusterCheck_Col*GetPitchU(),Nbins_ClusterCheck_Row*GetPitchV());
	
	double InCol_ppp  = -1.0;
	double OutCol_ppp = -1.0;
	double InRow_ppp  = -1.0;
	double OutRow_ppp = -1.0;
	ComputePixelPositionColRow_FromUV(InputUpos, InputVpos,InCol_ppp, InRow_ppp);
	ComputePixelPositionColRow_FromUV(OutputUpos,OutputVpos,OutCol_ppp,OutRow_ppp);
	
	InPostion[NhitsToDraw[ithr]][ithr] = new TEllipse(InCol_ppp,InRow_ppp,rell_ppp,rell_ppp);
	InPostion[NhitsToDraw[ithr]][ithr]->SetLineColor(1);
	InPostion[NhitsToDraw[ithr]][ithr]->SetLineWidth(1);
	InPostion[NhitsToDraw[ithr]][ithr]->SetFillColor(1);
	
	OutPostion[NhitsToDraw[ithr]][ithr] = new TEllipse(OutCol_ppp,OutRow_ppp,rell_ppp,rell_ppp);
	OutPostion[NhitsToDraw[ithr]][ithr]->SetLineColor(1);
	OutPostion[NhitsToDraw[ithr]][ithr]->SetLineWidth(1);
	OutPostion[NhitsToDraw[ithr]][ithr]->SetFillColor(6);
	
	InPostionUV[NhitsToDraw[ithr]][ithr] = new TEllipse(InputUpos, InputVpos,rell_pppUV,rell_pppUV);
	InPostionUV[NhitsToDraw[ithr]][ithr]->SetLineColor(1);
	InPostionUV[NhitsToDraw[ithr]][ithr]->SetLineWidth(1);
	InPostionUV[NhitsToDraw[ithr]][ithr]->SetFillColor(1);
	
	OutPostionUV[NhitsToDraw[ithr]][ithr] = new TEllipse(OutputUpos, OutputVpos,rell_pppUV,rell_pppUV);
	OutPostionUV[NhitsToDraw[ithr]][ithr]->SetLineColor(1);
	OutPostionUV[NhitsToDraw[ithr]][ithr]->SetLineWidth(1);
	OutPostionUV[NhitsToDraw[ithr]][ithr]->SetFillColor(6);
	
	NPixelsCheck[NhitsToDraw[ithr]][ithr] = 0;
	for(Int_t ipix = 0; ipix < int(fPixelMap.size()); ipix++){
	  double Decimals = 1.0e+3;
          Int_t col,row;
          GetXYPixelNumber(col,row,fPixelMap[ipix]);
	  Double_t u,v;
	  ComputePixelPositionUV_FromColRow(col,row,u,v);
	  
	  if(ipix+1 > MaxPixelInCluster) {
	    cout << "WARNING: for Threshold = " << ytitle << " (thr_idx = " << ithr+1 << "), evt = " << Nev+1 << ", the index of pixel in cluster " << ipix+1
	         << " goes beyong allowed limit " << MaxPixelInCluster << ". Cutting this cluster with size " << fPixelMap.size()
		 << endl;
	    continue;
	  }
 	  
	  PixelPosition[NhitsToDraw[ithr]][ithr][NPixelsCheck[NhitsToDraw[ithr]][ithr]]      = TVector2(u,v);
	  PixelAnalogCharge[NhitsToDraw[ithr]][ithr][NPixelsCheck[NhitsToDraw[ithr]][ithr]]  = RoundOff(Decimals*fAnalogChargeMap[ipix]*1.0e-3)/Decimals;
	  PixelDigitalCharge[NhitsToDraw[ithr]][ithr][NPixelsCheck[NhitsToDraw[ithr]][ithr]] = fDigitalChargeMap[ipix];
	  NPixelsCheck[NhitsToDraw[ithr]][ithr]++;
	  
	  int idx_col = h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetXaxis()->FindBin(col);
	  int idx_row = h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->GetYaxis()->FindBin(row);
	  h_ClusterShape_Analog[NhitsToDraw[ithr]][ithr]->SetBinContent(idx_col,idx_row, RoundOff(Decimals*fAnalogChargeMap[ipix]*1.0e-3)/Decimals);
	  
	  if(fDigitalChargeMap[ipix] < 1) continue;
	  h_ClusterShape_Digital[NhitsToDraw[ithr]][ithr]->SetBinContent(idx_col,idx_row,1.0);
	}
	
	NhitsToDraw[ithr]++;
      }
      
      if(verbose) cout << "Digital multiplicity = " << totalMult << endl;
      h_SeedCharge[ithr]->Fill(Seed_Charge*1.0e-3);
      if(Seed_Charge > Discri_Threshold_electrons) IsDetected_analog = true;
      
      if(TMath::Abs(SumOfWeights_analog) > 1.0e-8 && IsDetected_analog) {
	Position_analog_U /= SumOfWeights_analog;
	Position_analog_V /= SumOfWeights_analog;
	
	if(verbose) {
	  cout << " Delta Analog Position (u,v) = (" << AveUpos - Position_analog_U << "," << AveVpos - Position_analog_V << ") um" << endl;
	}
	
	h_DeltaAnalogPositionU[ithr]->Fill(AveUpos - Position_analog_U);
	h_DeltaAnalogPositionV[ithr]->Fill(AveVpos - Position_analog_V);
	h_DeltaAnalogPositionUV[ithr]->Fill(AveUpos - Position_analog_U,AveVpos - Position_analog_V);
      }
      if(totalMult > 0) {
	hmultiplicity[ithr]->Fill(totalMult);
	Ndet_evts++;
	int mult_idx = -999;
        for(int imult=0;imult<NMults;imult++) {
	  if(imult+1 < NMults) {
	    if(totalMult == imult+1) {
	      mult_idx = imult;
	      break;
	    }
	  }
	  else mult_idx = NMults-1;
        }
	
	Position_digital_U /= totalMult;
	Position_digital_V /= totalMult;
	
	if(verbose) {
	  cout << " Digital Position (u,v) = (" << AveUpos - Position_digital_U << "," << AveVpos - Position_digital_V << ") um" << endl;
	}
	
	h_DeltaDigitalPositionU[ithr]->Fill(AveUpos - Position_digital_U);
	h_DeltaDigitalPositionV[ithr]->Fill(AveVpos - Position_digital_V);
	h_DeltaDigitalPositionUV[ithr]->Fill(AveUpos - Position_digital_U,AveVpos - Position_digital_V);
	
	h_DeltaDigitalPositionU_vs_Mult[ithr][mult_idx]->Fill(AveUpos - Position_digital_U);
	h_DeltaDigitalPositionV_vs_Mult[ithr][mult_idx]->Fill(AveVpos - Position_digital_V);
	
	h_TrkDistToDiode_Detected[ithr]->Fill(trk_dist_to_diode);
	h_TrkDistToDiode_Detected_vs_Mult[ithr][mult_idx]->Fill(trk_dist_to_diode);
	
	h_TrkDistToDiode2D_Detected[ithr]->Fill(ShiftedPosition.X(),ShiftedPosition.Y());
	h_TrkDistToDiode2D_Detected_vs_Mult[ithr][mult_idx]->Fill(ShiftedPosition.X(),ShiftedPosition.Y());
      }
      else {
	h_TrkDistToDiode_NonDetected[ithr]->Fill(trk_dist_to_diode);
	h_TrkDistToDiode2D_NonDetected[ithr]->Fill(ShiftedPosition.X(),ShiftedPosition.Y());
      }

      if(verbose) {
        cout << "Event " << Nev+1 << " begin:"<< endl;
        cout << endl;
      }

    } // End of loops over events

    int MinNevts = 25;
    double efficiency,err_efficiency;
    double sigma_analog_U,err_sigma_analog_U;
    double sigma_analog_V,err_sigma_analog_V;
    double sigma_digital_U,err_sigma_digital_U;
    double sigma_digital_V,err_sigma_digital_V;
    if(Nevents > MinNevts) {
      efficiency     = (double)Ndet_evts/Nevents;
      err_efficiency = sqrt(efficiency*(1.0 - efficiency)/Nevents);
      
      efficiency     *= 100.0;
      err_efficiency *= 100.0;
      
      //cout << "efficiency = " << efficiency << " +/- " << err_efficiency << " %" << endl;
    }
    else {
      efficiency     = -1.0;
      err_efficiency = 1.0e-10;
    }
    h_effic_vs_Thr->SetBinContent(ithr+1,efficiency);
    h_effic_vs_Thr->SetBinError(ithr+1,err_efficiency);
    
    if(h_DeltaAnalogPositionU[ithr]->GetEntries() > MinNevts) {
      sigma_analog_U     = h_DeltaAnalogPositionU[ithr]->GetRMS();
      err_sigma_analog_U = h_DeltaAnalogPositionU[ithr]->GetRMSError();
    }
    else {
      sigma_analog_U     = -1.0;
      err_sigma_analog_U =  1.0e-10;
    }
    h_AnalogResolutionU_vs_Thr->SetBinContent(ithr+1,sigma_analog_U);
    h_AnalogResolutionU_vs_Thr->SetBinError(ithr+1,err_sigma_analog_U);
    
    if(h_DeltaAnalogPositionV[ithr]->GetEntries() > MinNevts) {
      sigma_analog_V     = h_DeltaAnalogPositionV[ithr]->GetRMS();
      err_sigma_analog_V = h_DeltaAnalogPositionV[ithr]->GetRMSError();
    }
    else {
      sigma_analog_V     = -1.0;
      err_sigma_analog_V =  1.0e-10;
    }
    h_AnalogResolutionV_vs_Thr->SetBinContent(ithr+1,sigma_analog_V);
    h_AnalogResolutionV_vs_Thr->SetBinError(ithr+1,err_sigma_analog_V);
    
    if(h_DeltaDigitalPositionU[ithr]->GetEntries() > MinNevts) {
      sigma_digital_U     = h_DeltaDigitalPositionU[ithr]->GetRMS();
      err_sigma_digital_U = h_DeltaDigitalPositionU[ithr]->GetRMSError();
    }
    else {
      sigma_digital_U     = -1.0;
      err_sigma_digital_U =  1.0e-10;
    }
    h_DigitalResolutionU_vs_Thr->SetBinContent(ithr+1,sigma_digital_U);
    h_DigitalResolutionU_vs_Thr->SetBinError(ithr+1,err_sigma_digital_U);
    
    if(h_DeltaDigitalPositionV[ithr]->GetEntries() > MinNevts) {
      sigma_digital_V     = h_DeltaDigitalPositionV[ithr]->GetRMS();
      err_sigma_digital_V = h_DeltaDigitalPositionV[ithr]->GetRMSError();
    }
    else {
      sigma_digital_V     = -1.0;
      err_sigma_digital_V =  1.0e-10;
    }
    h_DigitalResolutionV_vs_Thr->SetBinContent(ithr+1,sigma_digital_V);
    h_DigitalResolutionV_vs_Thr->SetBinError(ithr+1,err_sigma_digital_V);

  } // End loop over thresholds
  

  double Maximum,Minimum,porcent;

  TString EPSName;
  EPSName = TString(output_file) + TString(".eps");
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetRightMargin(0.15);
  
  c1->Print(EPSNameO.Data());
  
  for(int ithr=0;ithr<N_Thr_steps;ithr++) {
    TString FuncName;
    
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    hmultiplicity[ithr]->Draw();
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_SeedCharge[ithr]->Draw();
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_DeltaAnalogPositionUV[ithr]->Draw("colz");
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    
    FuncName = TString("Fit_") + TString(h_DeltaAnalogPositionV[ithr]->GetName());
    f_DeltaAnalogPositionV[ithr] = new TF1(HistName.Data(),FitModel.Data(),
					   -Fit_factor*h_DeltaAnalogPositionV[ithr]->GetRMS(),
					   +Fit_factor*h_DeltaAnalogPositionV[ithr]->GetRMS());
    f_DeltaAnalogPositionV[ithr]->SetLineColor(kRed);
    f_DeltaAnalogPositionV[ithr]->SetLineWidth(2);
    f_DeltaAnalogPositionV[ithr]->SetParameter(0,h_DeltaAnalogPositionV[ithr]->GetMean());
    f_DeltaAnalogPositionV[ithr]->SetParName(0,"#mu (#mum)");
    f_DeltaAnalogPositionV[ithr]->SetParameter(1,h_DeltaAnalogPositionV[ithr]->GetRMS());
    f_DeltaAnalogPositionV[ithr]->SetParName(1,"#sigma (#mum)");
    f_DeltaAnalogPositionV[ithr]->SetParameter(2,h_DeltaAnalogPositionV[ithr]->GetEntries());
    f_DeltaAnalogPositionV[ithr]->SetParName(2,"N_{evts}");
    
    h_DeltaAnalogPositionV[ithr]->Sumw2();
    h_DeltaAnalogPositionV[ithr]->Fit(f_DeltaAnalogPositionV[ithr],"R");
    h_DeltaAnalogPositionV[ithr]->Draw();
    c1->cd(3);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    
    FuncName = TString("Fit_") + TString(h_DeltaAnalogPositionU[ithr]->GetName());
    f_DeltaAnalogPositionU[ithr] = new TF1(HistName.Data(),FitModel.Data(),
					   -Fit_factor*h_DeltaAnalogPositionU[ithr]->GetRMS(),
					   +Fit_factor*h_DeltaAnalogPositionU[ithr]->GetRMS());
    f_DeltaAnalogPositionU[ithr]->SetLineColor(kRed);
    f_DeltaAnalogPositionU[ithr]->SetLineWidth(2);
    f_DeltaAnalogPositionU[ithr]->SetParameter(0,h_DeltaAnalogPositionU[ithr]->GetMean());
    f_DeltaAnalogPositionU[ithr]->SetParName(0,"#mu (#mum)");
    f_DeltaAnalogPositionU[ithr]->SetParameter(1,h_DeltaAnalogPositionU[ithr]->GetRMS());
    f_DeltaAnalogPositionU[ithr]->SetParName(1,"#sigma (#mum)");
    f_DeltaAnalogPositionU[ithr]->SetParameter(2,h_DeltaAnalogPositionU[ithr]->GetEntries());
    f_DeltaAnalogPositionU[ithr]->SetParName(2,"N_{evts}");
    
    h_DeltaAnalogPositionU[ithr]->Sumw2();
    h_DeltaAnalogPositionU[ithr]->Fit(f_DeltaAnalogPositionU[ithr],"R");
    h_DeltaAnalogPositionU[ithr]->Draw();
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_DeltaDigitalPositionUV[ithr]->Draw("colz");
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    
    FuncName = TString("Fit_") + TString(h_DeltaDigitalPositionV[ithr]->GetName());
    f_DeltaDigitalPositionV[ithr] = new TF1(HistName.Data(),FitModel.Data(),
					    -Fit_factor*h_DeltaDigitalPositionV[ithr]->GetRMS(),
					    +Fit_factor*h_DeltaDigitalPositionV[ithr]->GetRMS());
    f_DeltaDigitalPositionV[ithr]->SetLineColor(kRed);
    f_DeltaDigitalPositionV[ithr]->SetLineWidth(2);
    f_DeltaDigitalPositionV[ithr]->SetParameter(0,h_DeltaDigitalPositionV[ithr]->GetMean());
    f_DeltaDigitalPositionV[ithr]->SetParName(0,"#mu (#mum)");
    f_DeltaDigitalPositionV[ithr]->SetParameter(1,h_DeltaDigitalPositionV[ithr]->GetRMS());
    f_DeltaDigitalPositionV[ithr]->SetParName(1,"#sigma (#mum)");
    f_DeltaDigitalPositionV[ithr]->SetParameter(2,h_DeltaDigitalPositionV[ithr]->GetEntries());
    f_DeltaDigitalPositionV[ithr]->SetParName(2,"N_{evts}");
    
    h_DeltaDigitalPositionV[ithr]->Sumw2();
    h_DeltaDigitalPositionV[ithr]->Fit(f_DeltaDigitalPositionV[ithr],"R");
    h_DeltaDigitalPositionV[ithr]->Draw();
    c1->cd(3);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    
    FuncName = TString("Fit_") + TString(h_DeltaDigitalPositionU[ithr]->GetName());
    f_DeltaDigitalPositionU[ithr] = new TF1(HistName.Data(),FitModel.Data(),
					    -Fit_factor*h_DeltaDigitalPositionU[ithr]->GetRMS(),
					    +Fit_factor*h_DeltaDigitalPositionU[ithr]->GetRMS());
    f_DeltaDigitalPositionU[ithr]->SetLineColor(kRed);
    f_DeltaDigitalPositionU[ithr]->SetLineWidth(2);
    f_DeltaDigitalPositionU[ithr]->SetParameter(0,h_DeltaDigitalPositionU[ithr]->GetMean());
    f_DeltaDigitalPositionU[ithr]->SetParName(0,"#mu (#mum)");
    f_DeltaDigitalPositionU[ithr]->SetParameter(1,h_DeltaDigitalPositionU[ithr]->GetRMS());
    f_DeltaDigitalPositionU[ithr]->SetParName(1,"#sigma (#mum)");
    f_DeltaDigitalPositionU[ithr]->SetParameter(2,h_DeltaDigitalPositionU[ithr]->GetEntries());
    f_DeltaDigitalPositionU[ithr]->SetParName(2,"N_{evts}");
    
    h_DeltaDigitalPositionU[ithr]->Sumw2();
    h_DeltaDigitalPositionU[ithr]->Fit(f_DeltaDigitalPositionU[ithr],"R");
    h_DeltaDigitalPositionU[ithr]->Draw();
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,3);
    for(int imult=0;imult<NMults;imult++) {
      c1->cd(imult+1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
            
      FuncName = TString("Fit_") + TString(h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetName());
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult] = new TF1(HistName.Data(),FitModel.Data(),
							     -Fit_factor*h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetRMS(),
							     +Fit_factor*h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetRMS());
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetLineColor(kRed);
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetLineWidth(2);
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetParameter(0,h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetMean());
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetParName(0,"#mu (#mum)");
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetParameter(1,h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetRMS());
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetParName(1,"#sigma (#mum)");
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetParameter(2,h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetEntries());
      f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->SetParName(2,"N_{evts}");
      
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->Sumw2();
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->Fit(f_DeltaDigitalPositionU_vs_Mult[ithr][imult],"R");
      h_DeltaDigitalPositionU_vs_Mult[ithr][imult]->Draw();
      
      h_MeanDigitalDeltaU_vs_Mult[ithr]->SetBinContent(imult+1,f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetParameter(0));
      h_MeanDigitalDeltaU_vs_Mult[ithr]->SetBinError(imult+1,f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetParError(0));
      h_RMSDigitalDeltaU_vs_Mult[ithr]->SetBinContent(imult+1,f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetParameter(1));
      h_RMSDigitalDeltaU_vs_Mult[ithr]->SetBinError(imult+1,f_DeltaDigitalPositionU_vs_Mult[ithr][imult]->GetParError(1));
    }
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,3);
    for(int imult=0;imult<NMults;imult++) {
      c1->cd(imult+1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      
      FuncName = TString("Fit_") + TString(h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetName());
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult] = new TF1(HistName.Data(),FitModel.Data(),
							     -Fit_factor*h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetRMS(),
							     +Fit_factor*h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetRMS());
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetLineColor(kRed);
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetLineWidth(2);
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetParameter(0,h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetMean());
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetParName(0,"#mu (#mum)");
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetParameter(1,h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetRMS());
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetParName(1,"#sigma (#mum)");
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetParameter(2,h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetEntries());
      f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->SetParName(2,"N_{evts}");
      
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->Sumw2();
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->Fit(f_DeltaDigitalPositionV_vs_Mult[ithr][imult],"R");
      h_DeltaDigitalPositionV_vs_Mult[ithr][imult]->Draw();
      
      h_MeanDigitalDeltaV_vs_Mult[ithr]->SetBinContent(imult+1,f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetParameter(0));
      h_MeanDigitalDeltaV_vs_Mult[ithr]->SetBinError(imult+1,f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetParError(0));
      h_RMSDigitalDeltaV_vs_Mult[ithr]->SetBinContent(imult+1,f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetParameter(1));
      h_RMSDigitalDeltaV_vs_Mult[ithr]->SetBinError(imult+1,f_DeltaDigitalPositionV_vs_Mult[ithr][imult]->GetParError(1));
    }
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    Maximum = -1.0e+20;
    Minimum = +1.0e+20;
    porcent = 0.1;
    for(int i=0;i<h_MeanDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->GetNbins();i++) {
      double v,e;
      v = h_MeanDigitalDeltaU_vs_Mult[ithr]->GetBinContent(i+1);
      e = h_MeanDigitalDeltaU_vs_Mult[ithr]->GetBinError(i+1);
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
      
      v = h_MeanDigitalDeltaV_vs_Mult[ithr]->GetBinContent(i+1);
      e = h_MeanDigitalDeltaV_vs_Mult[ithr]->GetBinError(i+1);
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
    }
    h_MeanDigitalDeltaU_vs_Mult[ithr]->SetMinimum(Minimum - porcent*(Maximum - Minimum));
    h_MeanDigitalDeltaU_vs_Mult[ithr]->SetMaximum(Maximum + porcent*(Maximum - Minimum));
    h_MeanDigitalDeltaV_vs_Mult[ithr]->SetMinimum(Minimum - porcent*(Maximum - Minimum));
    h_MeanDigitalDeltaV_vs_Mult[ithr]->SetMaximum(Maximum + porcent*(Maximum - Minimum));
    h_MeanDigitalDeltaU_vs_Mult[ithr]->Draw();
    h_MeanDigitalDeltaV_vs_Mult[ithr]->Draw("same");
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    Maximum = -1.0e+20;
    Minimum = +1.0e+20;
    porcent = 0.1;
    for(int i=0;i<h_RMSDigitalDeltaU_vs_Mult[ithr]->GetXaxis()->GetNbins();i++) {
      double v,e;
      v = h_RMSDigitalDeltaU_vs_Mult[ithr]->GetBinContent(i+1);
      e = h_RMSDigitalDeltaU_vs_Mult[ithr]->GetBinError(i+1);
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
      
      v = h_RMSDigitalDeltaV_vs_Mult[ithr]->GetBinContent(i+1);
      e = h_RMSDigitalDeltaV_vs_Mult[ithr]->GetBinError(i+1);
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
    }
    h_RMSDigitalDeltaU_vs_Mult[ithr]->SetMinimum(Minimum - porcent*(Maximum - Minimum));
    h_RMSDigitalDeltaU_vs_Mult[ithr]->SetMaximum(Maximum + porcent*(Maximum - Minimum));
    h_RMSDigitalDeltaV_vs_Mult[ithr]->SetMinimum(Minimum - porcent*(Maximum - Minimum));
    h_RMSDigitalDeltaV_vs_Mult[ithr]->SetMaximum(Maximum + porcent*(Maximum - Minimum));
    h_RMSDigitalDeltaU_vs_Mult[ithr]->Draw();
    h_RMSDigitalDeltaV_vs_Mult[ithr]->Draw("same");
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_TrkDistToDiode2D_Detected[ithr]->Draw("colz");
    for(int i=0;i<4;i++) Diodes[i]->Draw("same");
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_TrkDistToDiode2D_NonDetected[ithr]->Draw("colz");
    for(int i=0;i<4;i++) Diodes[i]->Draw("same");
    c1->cd(3);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_TrkDistToDiode_Detected[ithr]->Draw();
    c1->cd(4);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_TrkDistToDiode_NonDetected[ithr]->Draw();
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,3);
    for(int imult=0;imult<NMults;imult++) {
      c1->cd(imult+1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      h_TrkDistToDiode2D_Detected_vs_Mult[ithr][imult]->Draw("colz");
      for(int i=0;i<4;i++) Diodes[i]->Draw("same");
    }
    c1->Print(EPSName.Data());
    
    c1->Clear();
    c1->Divide(2,3);
    for(int imult=0;imult<NMults;imult++) {
      c1->cd(imult+1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      h_TrkDistToDiode_Detected_vs_Mult[ithr][imult]->Draw("colz");
    }
    c1->Print(EPSName.Data());
    
    Maximum = -1.0e+20;
    Minimum = +1.0e+20;
    for(int i=0;i<h_TrkDistToDiode_Detected[ithr]->GetXaxis()->GetNbins();i++) {
      double Ndet   = h_TrkDistToDiode_Detected[ithr]->GetBinContent(i+1);
      double Nnodet = h_TrkDistToDiode_NonDetected[ithr]->GetBinContent(i+1);
      double Ntot   = Ndet + Nnodet;
      double effic,err_effic;
      if(Ntot > 25) {
	effic     = Ndet/Ntot;
	err_effic = sqrt(effic*(1.0 - effic)/Ntot);
	
	effic     *= 100.0;
	err_effic *= 100.0;
	
	if(Maximum < effic + err_effic) Maximum = effic + err_effic;
	if(Minimum > effic - err_effic) Minimum = effic - err_effic;
	
      }
      else {
	effic     = -1.0;
	err_effic = 1.0e-10;
      }
      h_Effic_vs_TrkDistToDiode[ithr]->SetBinContent(i+1,effic);
      h_Effic_vs_TrkDistToDiode[ithr]->SetBinError(i+1,err_effic);
    }
    porcent = 0.10;
    h_Effic_vs_TrkDistToDiode[ithr]->SetMinimum(Minimum - porcent*(Maximum - Minimum));
    h_Effic_vs_TrkDistToDiode[ithr]->SetMaximum(Maximum + porcent*(Maximum - Minimum));
    c1->Clear();
    h_Effic_vs_TrkDistToDiode[ithr]->Draw();
    c1->Print(EPSName.Data());
  }
  
  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.15);
  Maximum = -1.0e+20;
  Minimum = +1.0e+20;
  for(int i=0;i<h_effic_vs_Thr->GetXaxis()->GetNbins();i++) {
    double v,e;
    v = h_effic_vs_Thr->GetBinContent(i+1);
    e = h_effic_vs_Thr->GetBinError(i+1);
    if(v < 0.0) continue;
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.10;
  h_effic_vs_Thr->SetMinimum(Minimum - porcent*(Maximum - Minimum));
  h_effic_vs_Thr->SetMaximum(Maximum + porcent*(Maximum - Minimum));
  h_effic_vs_Thr->Draw("l");
  TLine* lh_effic = new TLine(h_effic_vs_Thr->GetXaxis()->GetXmin(),99.0,
                              h_effic_vs_Thr->GetXaxis()->GetXmax(),99.0);
  lh_effic->SetLineColor(1);
  lh_effic->SetLineWidth(1);
  lh_effic->SetLineStyle(2);
  lh_effic->Draw();
  Maximum = -1.0e+20;
  Minimum = +1.0e+20;
  for(int i=0;i<h_AnalogResolutionU_vs_Thr->GetXaxis()->GetNbins();i++) {
    double v,e;
    v = h_AnalogResolutionU_vs_Thr->GetBinContent(i+1);
    e = h_AnalogResolutionU_vs_Thr->GetBinError(i+1);
    if(v >= 0.0) {
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
    }
    v = h_AnalogResolutionV_vs_Thr->GetBinContent(i+1);
    e = h_AnalogResolutionV_vs_Thr->GetBinError(i+1);
    if(v >= 0.0) {
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
    }
    v = h_DigitalResolutionU_vs_Thr->GetBinContent(i+1);
    e = h_DigitalResolutionU_vs_Thr->GetBinError(i+1);
    if(v >= 0.0) {
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
    }
    v = h_DigitalResolutionV_vs_Thr->GetBinContent(i+1);
    e = h_DigitalResolutionV_vs_Thr->GetBinError(i+1);
    if(v >= 0.0) {
      if(Maximum < v+e) Maximum = v+e;
      if(Minimum > v-e) Minimum = v-e;
    }
    
  }
  porcent = 0.20;
  h_AnalogResolutionU_vs_Thr->SetMinimum(Minimum - porcent*(Maximum - Minimum));
  h_AnalogResolutionU_vs_Thr->SetMaximum(Maximum + porcent*(Maximum - Minimum));
  h_AnalogResolutionV_vs_Thr->SetMinimum(Minimum - porcent*(Maximum - Minimum));
  h_AnalogResolutionV_vs_Thr->SetMaximum(Maximum + porcent*(Maximum - Minimum));
  h_DigitalResolutionU_vs_Thr->SetMinimum(Minimum - porcent*(Maximum - Minimum));
  h_DigitalResolutionU_vs_Thr->SetMaximum(Maximum + porcent*(Maximum - Minimum));
  h_DigitalResolutionV_vs_Thr->SetMinimum(Minimum - porcent*(Maximum - Minimum));
  h_DigitalResolutionV_vs_Thr->SetMaximum(Maximum + porcent*(Maximum - Minimum));
  
  c1->cd(3);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.15);
  h_AnalogResolutionU_vs_Thr->Draw("l");
  h_AnalogResolutionV_vs_Thr->Draw("lsame");
  c1->cd(4);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.15);
  h_DigitalResolutionU_vs_Thr->Draw("l");
  h_DigitalResolutionV_vs_Thr->Draw("lsame");
  c1->cd(2);
  double margin = 0.2;
  TLegend* leg = new TLegend(margin,margin,1.0-margin,1.0-margin);
  leg->SetFillColor(10);
  leg->AddEntry(h_DigitalResolutionU_vs_Thr,"resolution U","l");
  leg->AddEntry(h_DigitalResolutionV_vs_Thr,"resolution V","l");
  leg->Draw();
  c1->Print(EPSName.Data());
  
  c1->Print(EPSNameC.Data());
  
  TString PDFName = TString(output_file) + TString(".pdf");
  TString command;
  command = TString("ps2pdf  ") + EPSName + TString("   ") + PDFName;
  gSystem->Exec(command.Data());
  command = TString("rm  ") + EPSName;
  gSystem->Exec(command.Data());
  
  
  EPSName  = TString(output_file) + TString("_ClusterCheck.eps");
  EPSNameO = EPSName + TString("[");
  EPSNameC = EPSName + TString("]");
  
  c1->Print(EPSNameO.Data());
  
  for(int ithr=0;ithr<N_Thr_steps;ithr++) {
    for(int ihit=0;ihit<NhitsToDraw[ithr];ihit++) {
      
      TLine lll(InPostion[ihit][ithr]->GetX1(), InPostion[ihit][ithr]->GetY1(),
	        OutPostion[ihit][ithr]->GetX1(),OutPostion[ihit][ithr]->GetY1());
      lll.SetLineColor(6);
      lll.SetLineWidth(2);
      lll.SetLineStyle(2);
      
      TLine lllUV(InPostionUV[ihit][ithr]->GetX1(), InPostionUV[ihit][ithr]->GetY1(),
	          OutPostionUV[ihit][ithr]->GetX1(),OutPostionUV[ihit][ithr]->GetY1());
      lllUV.SetLineColor(6);
      lllUV.SetLineWidth(2);
      lllUV.SetLineStyle(2);
      
      TLine vvv[h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetNbins()];
      TLine hhh[h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetNbins()];
      
      for(int i=0;i<h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetNbins()-1;i++) {
	double x = h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetBinCenter(i+1) + 0.5*h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetBinWidth(i+1);
	vvv[i] = TLine(x,h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetXmin(),
		       x,h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetXmax());
	vvv[i].SetLineColor(1);
	vvv[i].SetLineWidth(1);
	vvv[i].SetLineStyle(2);
      }
      for(int i=0;i<h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetNbins()-1;i++) {
	double x = h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetBinCenter(i+1) + 0.5*h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetBinWidth(i+1);
	hhh[i] = TLine(h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetXmin(),x,
		       h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetXmax(),x);
	hhh[i].SetLineColor(1);
	hhh[i].SetLineWidth(1);
	hhh[i].SetLineStyle(2);
      }
      
      TH2F* h_ClusterShape_AnalogUV[NPixelsCheck[ihit][ithr]];
      TH2F* h_ClusterShape_DigitalUV[NPixelsCheck[ihit][ithr]];
      
      double RangeClusterU[2];
      double RangeClusterV[2];
      RangeClusterU[0] = RangeClusterV[0] = +1.0e+20;
      RangeClusterU[1] = RangeClusterV[1] = -1.0e+20;
      for(int i=0;i<NPixelsCheck[ihit][ithr];i++) {
	double u = PixelPosition[ihit][ithr][i].X();
	double v = PixelPosition[ihit][ithr][i].Y();
	if(RangeClusterU[0] > u) RangeClusterU[0] = u;
	if(RangeClusterU[1] < u) RangeClusterU[1] = u;
	if(RangeClusterV[0] > v) RangeClusterV[0] = v;
	if(RangeClusterV[1] < v) RangeClusterV[1] = v;
	
	HistName = TString("h_ClusterShape_AnalogUV_Hit") + long(ihit+1) + TString("_Thr") + long(ithr+1) + TString("_pix") + long(i+1);
	h_ClusterShape_AnalogUV[i] = new TH2F(HistName.Data(),"",
					      1,u-0.5*GetPitchU(),u+0.5*GetPitchU(),
					      1,v-0.5*GetPitchV(),v+0.5*GetPitchV());
	h_ClusterShape_AnalogUV[i]->SetBinContent(1,1,PixelAnalogCharge[ihit][ithr][i]);
	
	HistName = TString("h_ClusterShape_DigitalUV_Hit") + long(ihit+1) + TString("_Thr") + long(ithr+1) + TString("_pix") + long(i+1);
	h_ClusterShape_DigitalUV[i] = new TH2F(HistName.Data(),"",
					      1,u-0.5*GetPitchU(),u+0.5*GetPitchU(),
					      1,v-0.5*GetPitchV(),v+0.5*GetPitchV());
	if(PixelDigitalCharge[ihit][ithr][i] > 0) h_ClusterShape_DigitalUV[i]->SetBinContent(1,1,1);
	else                                     h_ClusterShape_DigitalUV[i]->SetBinContent(1,1,-10);
	
      }
      RangeClusterU[0] -= 1.5*GetPitchU();
      RangeClusterU[1] += 1.5*GetPitchU();
      RangeClusterV[0] -= 1.5*GetPitchV();
      RangeClusterV[1] += 1.5*GetPitchV();
      
      //cout << "RU = (" << RangeClusterU[0] << "," << RangeClusterU[1] << ") um" << endl;
      //cout << "RV = (" << RangeClusterV[0] << "," << RangeClusterV[1] << ") um" << endl;
      
      HistName = TString("h_ref_Hit") + long(ihit+1) + TString("_Thr") + long(ithr+1);
      TH2F h_ref(HistName.Data(),"",100,RangeClusterU[0],RangeClusterU[1],100,RangeClusterV[0],RangeClusterV[1]);
      h_ref.SetXTitle("U (#mum)");
      h_ref.GetXaxis()->CenterTitle(true);
      h_ref.SetYTitle("V (#mum)");
      h_ref.GetYaxis()->CenterTitle(true);
      h_ref.SetZTitle("");
      h_ref.GetZaxis()->CenterTitle(true);
      h_ref.GetXaxis()->SetTitleSize(TheSize);
      h_ref.GetXaxis()->SetLabelSize(TheSize);
      h_ref.GetYaxis()->SetTitleSize(TheSize);
      h_ref.GetYaxis()->SetLabelSize(TheSize);
      h_ref.GetZaxis()->SetTitleSize(TheSize);
      h_ref.GetZaxis()->SetLabelSize(TheSize);
      h_ref.SetStats(false);
      h_ref.SetBinContent(1,1,-1.0e+10);
      
      //PixelPosition[NhitsToDraw[ithr]][ithr][NPixelsCheck[NhitsToDraw[ithr]][ithr]]      = TVector2(u,v);
      //PixelAnalogCharge[NhitsToDraw[ithr]][ithr][NPixelsCheck[NhitsToDraw[ithr]][ithr]]  = RoundOff(Decimals*fAnalogChargeMap[ipix]*1.0e-3)/Decimals;
      //PixelDigitalCharge[NhitsToDraw[ithr]][ithr][NPixelsCheck[NhitsToDraw[ithr]][ithr]] = fDigitalChargeMap[ipix];
      //NPixelsCheck[NhitsToDraw[ithr]][ithr]++;
      //InPostionUV[NhitsToDraw[ithr]][ithr]->SetFillColor(1);
      
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      h_ClusterShape_Analog[ihit][ithr]->SetMarkerSize(1.8);
      h_ClusterShape_Analog[ihit][ithr]->Draw("textcolz");
      lll.Draw();
      InPostion[ihit][ithr]->Draw("same");
      OutPostion[ihit][ithr]->Draw("same");
      for(int i=0;i<h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetNbins()-1;i++) vvv[i].Draw();
      for(int i=0;i<h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetNbins()-1;i++) hhh[i].Draw();
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      h_ClusterShape_Digital[ihit][ithr]->Draw("col");
      lll.Draw();
      InPostion[ihit][ithr]->Draw("same");
      OutPostion[ihit][ithr]->Draw("same");
      for(int i=0;i<h_ClusterShape_Analog[ihit][ithr]->GetXaxis()->GetNbins()-1;i++) vvv[i].Draw();
      for(int i=0;i<h_ClusterShape_Analog[ihit][ithr]->GetYaxis()->GetNbins()-1;i++) hhh[i].Draw();
      c1->cd(3);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      h_ref.SetMinimum(h_ClusterShape_Analog[ihit][ithr]->GetMinimum());
      h_ref.SetMaximum(h_ClusterShape_Analog[ihit][ithr]->GetMaximum());
      //cout << "Min = " << h_ref.GetMinimum() << ", Max = " << h_ref.GetMaximum() << endl;
      h_ref.SetTitle("Analog 2D cluster");
      h_ref.SetZTitle("Analog Charge (#times10^{3} electrons)");
      h_ref.DrawClone("colz");
      for(int i=0;i<NPixelsCheck[ihit][ithr];i++) {
	h_ClusterShape_AnalogUV[i]->SetMarkerSize(1.8);
	h_ClusterShape_AnalogUV[i]->SetMinimum(h_ref.GetMinimum());
	h_ClusterShape_AnalogUV[i]->SetMaximum(h_ref.GetMaximum());
	h_ClusterShape_AnalogUV[i]->Draw("textsamecol");
      }
      lllUV.Draw();
      InPostionUV[ihit][ithr]->Draw("same");
      OutPostionUV[ihit][ithr]->Draw("same");
      c1->cd(4);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.15);
      h_ref.SetMinimum(1.0e-3);
      h_ref.SetMaximum(1.0 + 1.0e-3);
      h_ref.SetTitle("Digital 2D cluster");
      h_ref.SetZTitle("");
      h_ref.DrawClone("col");
      for(int i=0;i<NPixelsCheck[ihit][ithr];i++) {
	h_ClusterShape_DigitalUV[i]->SetMarkerSize(1.8);
	h_ClusterShape_DigitalUV[i]->SetMinimum(h_ref.GetMinimum());
	h_ClusterShape_DigitalUV[i]->SetMaximum(h_ref.GetMaximum());
	h_ClusterShape_DigitalUV[i]->Draw("samecol");
      }
      lllUV.Draw();
      InPostionUV[ihit][ithr]->Draw("same");
      OutPostionUV[ihit][ithr]->Draw("same");
      c1->Print(EPSName.Data());
      
      for(int i=0;i<NPixelsCheck[ihit][ithr];i++) {
	delete h_ClusterShape_AnalogUV[i];
      }
      
    }
  }
  c1->Print(EPSNameC.Data());
  
  PDFName = TString(output_file) + TString("_ClusterCheck.pdf");
  command = TString("ps2pdf  ") + EPSName + TString("   ") + PDFName;
  gSystem->Exec(command.Data());
  command = TString("rm  ") + EPSName;
  gSystem->Exec(command.Data());
  
  delete c1;
  delete r3;

  return;
  
}
//______________________________________________________________________________
//

