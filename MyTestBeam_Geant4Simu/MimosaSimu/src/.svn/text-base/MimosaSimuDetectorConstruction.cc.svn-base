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
// $Id: MimosaSimuDetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file MimosaSimuDetectorConstruction.cc
/// \brief Implementation of the MimosaSimuDetectorConstruction class
 
#include "MimosaSimuDetectorConstruction.hh"
#include "MimosaSimuDetectorMessenger.hh"
#include "MimosaSimuTrackerSD.hh"
#include "MimosaSimuTrajectorySD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* MimosaSimuDetectorConstruction::fMagFieldMessenger = 0;

MimosaSimuDetectorConstruction::MimosaSimuDetectorConstruction(MimosaSimuSetup* TheSetup) : G4VUserDetectorConstruction(), 
                                                                                            fNumberOfSensors(0),
                                                                                            fNbOfChambers(0),
                                                                                            fLogicMimosaSensor(NULL),
                                                                                            fTargetMaterial(NULL), 
                                                                                            fChamberMaterial(NULL),
                                                                                            fSensorMaterial(NULL),
                                                                                            fStepLimit(NULL),
                                                                                            fCheckOverlaps(true),
                                                                                            fSetup(TheSetup)
{
  
  fNumberOfSensors             = fSetup->TrackerParameter.Planes;
  fSensorMaterial              = new G4Material*[fNumberOfSensors];
  fLogicMimosaSensor           = new G4LogicalVolume*[fNumberOfSensors];
  fPhysicMimosaSensor          = new G4VPhysicalVolume*[fNumberOfSensors];
  
  fLogicMimosaSensorMetal      = new G4LogicalVolume*[fNumberOfSensors];
  fPhysicMimosaSensorMetal     = new G4VPhysicalVolume*[fNumberOfSensors];
  fLogicMimosaSensorSubstrate  = new G4LogicalVolume*[fNumberOfSensors];
  fPhysicMimosaSensorSubstrate = new G4VPhysicalVolume*[fNumberOfSensors];
  
  fLogicSource                 = NULL;
  fPhysicSource                = NULL;
  
  fMessenger = new MimosaSimuDetectorMessenger(this);
  
  //Setting up the name of the output GDML file
  TString GDMLOutPutName = fSetup->GetConfigFileName() + TString(".gdml");
  TString command = TString("rm ") + GDMLOutPutName;
  gSystem->Exec(command.Data());
  fWriteFile = GDMLOutPutName.Data();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
MimosaSimuDetectorConstruction::~MimosaSimuDetectorConstruction()
{
  delete [] fLogicMimosaSensor;
  delete [] fLogicMimosaSensorMetal;
  delete [] fLogicMimosaSensorSubstrate;
  
  delete [] fPhysicMimosaSensor;
  delete [] fPhysicMimosaSensorMetal;
  delete [] fPhysicMimosaSensorSubstrate;
  
  delete [] fSensorMaterial; 
  
  if(fStepLimit) {
    delete fStepLimit;
    fStepLimit = NULL;
  }
  if(fMessenger) {
    delete fMessenger;
    fMessenger = NULL;
  }
  
  //cout << "MimosaSimuDetectorConstruction::~MimosaSimuDetectorConstruction:: All OK!!!" << endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* MimosaSimuDetectorConstruction::Construct()
{
  
  G4VPhysicalVolume* fWorldPhysVol;
  
  // Define materials
  DefineMaterials();

  // Define volumes
  fWorldPhysVol = DefineVolumes();
  
  //Write geometry to GDML file
  fParser.Write(fWriteFile,fWorldPhysVol);
  
  return fWorldPhysVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();
  
  air        = nistManager->FindOrBuildMaterial("G4_AIR");
  silicon    = nistManager->FindOrBuildMaterial("G4_Si");
  carbon     = nistManager->FindOrBuildMaterial("G4_C");
  beryllium  = nistManager->FindOrBuildMaterial("G4_Be");
  tape       = nistManager->FindOrBuildMaterial("G4_POLYVINYL_ACETATE");
  mylar      = nistManager->FindOrBuildMaterial("G4_MYLAR");
  vacuum     = nistManager->FindOrBuildMaterial("G4_Galactic");
  aluminium  = nistManager->FindOrBuildMaterial("G4_Al");
  copper     = nistManager->FindOrBuildMaterial("G4_Cu");
  kapton     = nistManager->FindOrBuildMaterial("G4_KAPTON");
  
  //PLUME Silicon-Carbide composition
  double Non_air_fraction = 4.0*perCent;
  double fraction_C       = 0.5;
  double density_Si  = silicon->GetDensity();
  double density_C   = carbon->GetDensity();
  double density_air = air->GetDensity();
  double Frac_Si     = (1.0 - fraction_C) * Non_air_fraction;
  double Frac_C      = fraction_C         * Non_air_fraction;
  double Frac_Air    = 1.0 - Frac_Si - Frac_C;
  double density_SiC = density_air*Frac_Air + density_Si*Frac_Si + density_C*Frac_C;
  SiC_foam   = new G4Material("SiliconCarbide",density_SiC,3);
  SiC_foam->AddMaterial(silicon, Frac_Si);
  SiC_foam->AddMaterial(carbon,  Frac_C);
  SiC_foam->AddMaterial(air,     Frac_Air);
  
  glue     = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass"); //Need to find a better model of epoxy glue
  
  //if(!vacuum) cout << "Didn't find vaccum material" << endl; 
  
  //Filling up the material map
  _material_name_list.clear();
  _material_map[TString("silicon")] = silicon;
  _material_name_list.push_back(TString("silicon"));
  _material_map[TString("tape")]    = tape;
  _material_name_list.push_back(TString("tape"));
  _material_map[TString("mylar")]    = mylar;
  _material_name_list.push_back(TString("mylar"));
  _material_map[TString("Vacuum")]    = vacuum;
  _material_name_list.push_back(TString("Vacuum"));
  _material_map[TString("air")]    = air;
  _material_name_list.push_back(TString("air"));
  _material_map[TString("aluminium")]    = aluminium;
  _material_name_list.push_back(TString("aluminium"));
  _material_map[TString("carbon")]       = carbon;
  _material_name_list.push_back(TString("carbon"));
  _material_map[TString("beryllium")]    = beryllium;
  _material_name_list.push_back(TString("beryllium"));
  _material_map[TString("copper")]       = copper;
  _material_name_list.push_back(TString("copper"));
  _material_map[TString("kapton")]       = kapton;
  _material_name_list.push_back(TString("kapton"));
  _material_map[TString("SiC")]          = SiC_foam;
  _material_name_list.push_back(TString("SiC"));
  _material_map[TString("glue")]          = glue;
  _material_name_list.push_back(TString("glue"));
  
  //G4double a, z,density;
  //silicon = new G4Material("silicon", z=14., a=  28.09*g/mole, density=2.3290*g/cm3);
  
  // Lead defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");

  // Xenon gas defined using NIST Manager
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Xe");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MimosaSimuDetectorConstruction::DefineVolumes()
{

  G4String Some_LV;
  
  // Sizes of the principal geometrical components (solids)
  
  //Source parameters
  G4ThreeVector SourcePosition = G4ThreeVector(0.0,0.0,0.0);
  G4ThreeVector SourceTilt     = G4ThreeVector(0.0,0.0,0.0);
  G4RotationMatrix *Source_rot    = NULL;
  G4RotationMatrix *Source_Invrot = NULL;
  G4double      SourceRadius   = 0.0;
  G4double      SourceHight    = 0.5*mm;
  G4String      SourceMaterial = TString("aluminium").Data();
  G4String      Source_Solid   = TString("Source_Solid").Data();
  G4String      Source_LV      = TString("Source_LV").Data();
  G4String      Source_PV      = TString("Source_PV").Data();
  
  // Mimosa Sensors 
  //const int MyNumberOfSensors(50);
  G4double fSensorPitchX[MyNumberOfSensors];
  G4double fSensorPitchY[MyNumberOfSensors];
  G4int    colNumberSensor[MyNumberOfSensors];
  G4int    rowNumberSensor[MyNumberOfSensors];
  G4double fSensorLengthX[MyNumberOfSensors];
  G4double fSensorLengthY[MyNumberOfSensors];
  G4int    ReadOut[MyNumberOfSensors];

  G4double fSensorThickness[MyNumberOfSensors];
  TString  fMaterial[MyNumberOfSensors];
  G4double fSensorMetalThickness[MyNumberOfSensors];
  G4double fSensorEpiThickness[MyNumberOfSensors];

  G4double SensorX[MyNumberOfSensors];          //mm
  G4double SensorY[MyNumberOfSensors];          //mm
  G4double SensorZ[MyNumberOfSensors];          //mm
  G4double SensorRotationX[MyNumberOfSensors];  //rad
  G4double SensorRotationY[MyNumberOfSensors];  //rad
  G4double SensorRotationZ[MyNumberOfSensors];  //rad
  G4ThreeVector SensorPosition[MyNumberOfSensors];
  G4RotationMatrix *rot[MyNumberOfSensors];
  G4RotationMatrix *Invrot[MyNumberOfSensors];

  G4double Rad_to_Deg = 180.0/TMath::Pi();
  
  G4double RX[2];
  G4double RY[2];
  G4double RZ[2];
  RX[0] = RY[0] = RZ[0] = +1.0e+20;
  RX[1] = RY[1] = RZ[1] = -1.0e+20;
  for(G4int sensorID = 0; sensorID < fSetup->TrackerParameter.Planes; sensorID++) { //Looping over the planes to build the
    //The sensor size:
    fSensorPitchX[sensorID]         = fSetup->GetPlanePar(sensorID+1).Pitch[0]*um;
    fSensorPitchY[sensorID]         = fSetup->GetPlanePar(sensorID+1).Pitch[1]*um;
    colNumberSensor[sensorID]       = fSetup->GetPlanePar(sensorID+1).Strips[0];
    rowNumberSensor[sensorID]       = fSetup->GetPlanePar(sensorID+1).Strips[1];
    fSensorLengthX[sensorID]        = colNumberSensor[sensorID]*fSensorPitchX[sensorID];
    fSensorLengthY[sensorID]        = rowNumberSensor[sensorID]*fSensorPitchY[sensorID];
    ReadOut[sensorID]               = fSetup->GetPlanePar(sensorID+1).Readout;
    
    //The sensor thickness and material
    fSensorThickness[sensorID]      = fSetup->GetPlanePar(sensorID+1).PlaneThickness*um;
    fSensorMetalThickness[sensorID] = fSetup->GetPlanePar(sensorID+1).PlaneMetalThickness;
    fSensorEpiThickness[sensorID]   = fSetup->GetPlanePar(sensorID+1).PlaneEpiThickness;
    fMaterial[sensorID]             = fSetup->GetPlanePar(sensorID+1).PlaneMaterial;

    //The sensor position
    SensorX[sensorID]               = fSetup->GetPlanePar(sensorID+1).Position[0]*um;
    SensorY[sensorID]               = fSetup->GetPlanePar(sensorID+1).Position[1]*um;
    SensorZ[sensorID]               = fSetup->GetPlanePar(sensorID+1).Position[2]*um;
    SensorRotationZ[sensorID]       = fSetup->GetPlanePar(sensorID+1).Tilt[0]*Rad_to_Deg*deg;
    SensorRotationY[sensorID]       = fSetup->GetPlanePar(sensorID+1).Tilt[1]*Rad_to_Deg*deg;
    SensorRotationX[sensorID]       = fSetup->GetPlanePar(sensorID+1).Tilt[2]*Rad_to_Deg*deg;
    
    fSensorMaterial[sensorID]       = silicon;
    
    SensorPosition[sensorID] = G4ThreeVector(SensorX[sensorID],SensorY[sensorID],SensorZ[sensorID]);
    
    rot[sensorID]    = new G4RotationMatrix();
    Invrot[sensorID] = new G4RotationMatrix();

    *rot[sensorID]    = G4RotationMatrix();
    *Invrot[sensorID] = G4RotationMatrix();
    
    G4RotationMatrix* rotX = new G4RotationMatrix();
    G4RotationMatrix* rotY = new G4RotationMatrix();
    G4RotationMatrix* rotZ = new G4RotationMatrix();
    *rotX = G4RotationMatrix();
    *rotY = G4RotationMatrix();
    *rotZ = G4RotationMatrix();
    rotX->rotateX(-SensorRotationX[sensorID]);
    rotY->rotateY(-SensorRotationY[sensorID]);
    rotZ->rotateZ(-SensorRotationZ[sensorID]);
    
    //*rot[sensorID]  = *rotX;
    //*rot[sensorID] *= *rotY;
    //*rot[sensorID] *= *rotZ;
    *Invrot[sensorID]  = *rotX;
    *Invrot[sensorID] *= *rotY;
    *Invrot[sensorID] *= *rotZ;
    
    //SetRotationMatrix(-SensorRotationX[sensorID],-SensorRotationY[sensorID],-SensorRotationZ[sensorID],rot[sensorID]);
    //rot[sensorID]->rotateX(-SensorRotationX[sensorID]);
    //rot[sensorID]->rotateY(-SensorRotationY[sensorID]);
    //rot[sensorID]->rotateZ(-SensorRotationZ[sensorID]);
    //*Invrot[sensorID] = rot[sensorID]->inverse();
    *rot[sensorID] = Invrot[sensorID]->inverse();

    //Building the dimensions of the world volume
    //Check that coordinates of the of all the corners of the Plane
    //Use the maximum/minimum X,Y,Z positions to build a box for the world volume
    G4RotationMatrix* temp_rot = rot[sensorID];
    //G4RotationMatrix* temp_rot = Invrot[sensorID];
    G4ThreeVector delta;
    //Lower-left corner
    delta = G4ThreeVector(-0.5*fSensorLengthX[sensorID],
			  -0.5*fSensorLengthY[sensorID],
			  +0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    delta = G4ThreeVector(-0.5*fSensorLengthX[sensorID],
			  -0.5*fSensorLengthY[sensorID],
			  -0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    
    //Upper-left corner
    delta = G4ThreeVector(-0.5*fSensorLengthX[sensorID],
			  +0.5*fSensorLengthY[sensorID],
			  +0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    delta = G4ThreeVector(-0.5*fSensorLengthX[sensorID],
			  +0.5*fSensorLengthY[sensorID],
			  -0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    
    //Lower-right corner
    delta = G4ThreeVector(+0.5*fSensorLengthX[sensorID],
			  -0.5*fSensorLengthY[sensorID],
			  +0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    delta = G4ThreeVector(+0.5*fSensorLengthX[sensorID],
			  -0.5*fSensorLengthY[sensorID],
			  -0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    
    //Upper-right corner
    delta = G4ThreeVector(+0.5*fSensorLengthX[sensorID],
			  +0.5*fSensorLengthY[sensorID],
			  +0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
    delta = G4ThreeVector(+0.5*fSensorLengthX[sensorID],
			  +0.5*fSensorLengthY[sensorID],
			  -0.5*fSensorThickness[sensorID]);
    delta.transform(*temp_rot);
    delta += SensorPosition[sensorID];
    if(RX[0] > delta(0)) RX[0] = delta(0);
    if(RX[1] < delta(0)) RX[1] = delta(0);
    if(RY[0] > delta(1)) RY[0] = delta(1);
    if(RY[1] < delta(1)) RY[1] = delta(1);
    if(RZ[0] > delta(2)) RZ[0] = delta(2);
    if(RZ[1] < delta(2)) RZ[1] = delta(2);
  }
  
  cout << "Range X = (" << RX[0]/mm << "," << RX[1]/mm << ") cm" << endl;
  cout << "Range Y = (" << RY[0]/mm << "," << RY[1]/mm << ") cm" << endl;
  cout << "Range Z = (" << RZ[0]/mm << "," << RZ[1]/mm << ") cm" << endl;
  
  if(fSetup->TrackerParameter.ExpSetup == TString("Beam-Test")) {
    if(!(fSetup->TrackerParameter.BeamOrigin == G4ThreeVector(0.0,0.0,-1.0e+6) || fSetup->TrackerParameter.BeamOriginSpread(0) < 0.0 || fSetup->TrackerParameter.BeamOriginSpread(1) < 0.0)) {
      //Checking the beam particles origin position
      G4ThreeVector delta  = fSetup->TrackerParameter.BeamOrigin;
      delta *= um;
      //cout << "Beam origin position = (" << delta(0)/cm << "," << delta(1)/cm << delta(2)/cm << ") cm" << endl;
      if(RX[0] > delta(0)) RX[0] = delta(0);
      if(RX[1] < delta(0)) RX[1] = delta(0);
      if(RY[0] > delta(1)) RY[0] = delta(1);
      if(RY[1] < delta(1)) RY[1] = delta(1);
      if(RZ[0] > delta(2)) RZ[0] = delta(2);
      if(RZ[1] < delta(2)) RZ[1] = delta(2);
    }
  }
  else if(fSetup->CheckIfSourceExists(fSetup->TrackerParameter.ExpSetup)) {
  //else if(false) {

    SourcePosition = G4ThreeVector(fSetup->TrackerParameter.SourcePosition(0)*um,
                                   fSetup->TrackerParameter.SourcePosition(1)*um,
				   fSetup->TrackerParameter.SourcePosition(2)*um);
    SourceTilt     = G4ThreeVector(fSetup->TrackerParameter.SourceTilt(2),
                                   fSetup->TrackerParameter.SourceTilt(1),
				   fSetup->TrackerParameter.SourceTilt(0));
    Source_rot    = new G4RotationMatrix();
    Source_Invrot = new G4RotationMatrix();

    *Source_rot = G4RotationMatrix();
    SetRotationMatrix(SourceTilt(0),SourceTilt(1),SourceTilt(2),Source_rot);
    *Source_Invrot = Source_rot->inverse();
    
    SourceRadius  = fSetup->TrackerParameter.SourceRadius*um;
    
    cout << endl;
    cout << "Building source volume:" << endl;
    cout << " - Cylinder, radius = " << SourceRadius/mm << " mm, and hight = " << SourceHight/mm << " mm" << endl;
    cout << " - Source position = (" << SourcePosition(0)/mm << "," << SourcePosition(1)/mm << "," <<  SourcePosition(2)/mm << ") mm"  << endl;
    cout << " - Source tilt     = (" << SourceTilt(2)/deg    << "," << SourceTilt(1)/deg    << "," <<  SourceTilt(0)/deg    << ") deg" << endl;
    cout << endl;

    G4RotationMatrix* temp_rot = Source_rot;
    G4int Nangle_scan_points = 10;
    G4ThreeVector ScanPosition;
    for(int ipoint=0;ipoint<Nangle_scan_points;ipoint++) {
      G4double phi_scan = (ipoint+0.5)*2.0*TMath::Pi()/Nangle_scan_points;
      
      ScanPosition = G4ThreeVector(SourceRadius*cos(phi_scan),
	                           SourceRadius*sin(phi_scan),
				   -0.5*SourceHight);
      ScanPosition.transform(*temp_rot);
      ScanPosition += SourcePosition;
      if(RX[0] > ScanPosition(0)) RX[0] = ScanPosition(0);
      if(RX[1] < ScanPosition(0)) RX[1] = ScanPosition(0);
      if(RY[0] > ScanPosition(1)) RY[0] = ScanPosition(1);
      if(RY[1] < ScanPosition(1)) RY[1] = ScanPosition(1);
      if(RZ[0] > ScanPosition(2)) RZ[0] = ScanPosition(2);
      if(RZ[1] < ScanPosition(2)) RZ[1] = ScanPosition(2);
      
      ScanPosition = G4ThreeVector(SourceRadius*cos(phi_scan),
	                           SourceRadius*sin(phi_scan),
				   +0.5*SourceHight);
      ScanPosition.transform(*temp_rot);
      ScanPosition += SourcePosition;
      if(RX[0] > ScanPosition(0)) RX[0] = ScanPosition(0);
      if(RX[1] < ScanPosition(0)) RX[1] = ScanPosition(0);
      if(RY[0] > ScanPosition(1)) RY[0] = ScanPosition(1);
      if(RY[1] < ScanPosition(1)) RY[1] = ScanPosition(1);
      if(RZ[0] > ScanPosition(2)) RZ[0] = ScanPosition(2);
      if(RZ[1] < ScanPosition(2)) RZ[1] = ScanPosition(2);
    }

  }
  else {}
  cout << "Range X = (" << RX[0]/mm << "," << RX[1]/mm << ") cm" << endl;
  cout << "Range Y = (" << RY[0]/mm << "," << RY[1]/mm << ") cm" << endl;
  cout << "Range Z = (" << RZ[0]/mm << "," << RZ[1]/mm << ") cm" << endl;
  
  G4double porcent = 20.0/100.0;
  G4double delta_range = 0.0;
  delta_range = RX[1] - RX[0];
  RX[0] -= delta_range*porcent;
  RX[1] += delta_range*porcent;
  delta_range = RY[1] - RY[0];
  RY[0] -= delta_range*porcent;
  RY[1] += delta_range*porcent;
  delta_range = RZ[1] - RZ[0];
  RZ[0] -= delta_range*porcent;
  RZ[1] += delta_range*porcent;
  
  G4double Xmax = TMath::Max(TMath::Abs(RX[0]),TMath::Abs(RX[1]));
  G4double Ymax = TMath::Max(TMath::Abs(RY[0]),TMath::Abs(RY[1]));
  G4double Zmax = TMath::Max(TMath::Abs(RZ[0]),TMath::Abs(RZ[1]));
  //Xmax = TMath::Max(RX[1],RY[1]);
  //Ymax = TMath::Max(RX[1],RY[1]);
  Xmax *= 2.0;
  Ymax *= 2.0;
  Zmax *= 2.0;
  
  if(TMath::Abs(Zmax/TMath::Max(Xmax,Ymax)) < 1.0) Zmax = TMath::Max(Xmax,Ymax);
  
  
  //TEST values
  //Xmax = 2.0*cm;
  //Ymax = 2.0*cm;
  //Zmax = 2.0*cm;
  
  cout << "Xmax = " << Xmax/cm << " cm" << endl;
  cout << "Ymax = " << Ymax/cm << " cm" << endl;
  cout << "Zmax = " << Zmax/cm << " cm" << endl;
  
  fWorldLengthX = TMath::Max(Xmax,Ymax);
  fWorldLengthY = TMath::Max(Xmax,Ymax);
  fWorldLengthZ = Zmax;
  cout << "fWorldLengthX/2 = " << 0.5*fWorldLengthX/cm << " cm" << endl;
  cout << "fWorldLengthY/2 = " << 0.5*fWorldLengthY/cm << " cm" << endl;
  cout << "fWorldLengthZ/2 = " << 0.5*fWorldLengthZ/cm << " cm" << endl;
  G4ThreeVector fPositionWorld(0,0,0);
  G4double MaxworldLength = TMath::Max(fWorldLengthX,fWorldLengthY);
  MaxworldLength = TMath::Max(MaxworldLength,fWorldLengthZ);
  
  // World
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(MaxworldLength);
  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;
	 
  G4Box* worldS = new G4Box("world",                                    //its name
                            fWorldLengthX/2,
			    fWorldLengthY/2,
			    fWorldLengthZ/2); //its size
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,   //its solid
						 _material_map[fSetup->TrackerParameter.MediumMaterial], //its material
						 "World"); //its name
						 
  Some_LV = "World";
  _ListOfAllLogicalVolume_Names.push_back(Some_LV);
  _ListOfAllLogicalVolumes.push_back(worldLV);
  G4double World_maxStep = TMath::Min(fWorldLengthX,fWorldLengthY);
  World_maxStep  = TMath::Min(World_maxStep,fWorldLengthY);
  if(fSetup->TrackerParameter.BFieldMagnitude < 1.0e-3) World_maxStep *= 1.0e-1;
  else                                                  World_maxStep *= 1.0e-2;
  worldLV->SetUserLimits(new G4UserLimits(World_maxStep)); 
  
  //  Must place the World Physical volume unrotated at (0,0,0)
  G4VPhysicalVolume* worldPV = new G4PVPlacement(0,               // no rotation
						 fPositionWorld,  // at (0,0,0)
                                                 worldLV,         // its logical volume
                                                 "World",         // its name
                                                 0,               // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps 
  
  // Definitions of Solids, Logical Volumes, Physical Volumes

  // Visualization attributes
  G4VisAttributes* boxVisAtt              = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* MimosaVisAtt           = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* MimosaMetalVisAtt      = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  G4VisAttributes* MimosaSubstrateVisAtt  = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* InnertVisAtt           = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* TheVisAttSource        = new G4VisAttributes(G4Colour(1.0,0.0,1.0));

  worldLV->SetVisAttributes(boxVisAtt);

  for(G4int sensorID = 0; sensorID < fSetup->TrackerParameter.Planes; sensorID++) { //Looping over the planes to build the
    cout << "Sensor " << sensorID+1 << ": location = (" << SensorX[sensorID]/cm << "," << SensorY[sensorID]/cm << "," << SensorZ[sensorID]/cm << ") cm"
         << "; SizeU/V/W = (" << fSensorLengthX[sensorID]/cm << " cm," << fSensorLengthY[sensorID]/cm << " cm," << fSensorThickness[sensorID]/um << " um)"
	 << "; rotation X/Y/Z = (" << SensorRotationX[sensorID]/deg << "," << SensorRotationY[sensorID]/deg << "," << SensorRotationZ[sensorID]/deg << ") deg"
         << endl;

    //notActivePartTop[sensorID]    = 0.0; // um
    //notActivePartBottom[sensorID] = 0.0; // um
	 
    bool WithEpiMetalAndSubstrate = false;
    if(1.0 - fSensorEpiThickness[sensorID] > 1.0e-3) WithEpiMetalAndSubstrate = true;

    G4String Solid_Name;
    G4String LV_Name;
    G4String PV_Name;
    
    G4String Solid_Name_metal;
    G4String LV_Name_metal;
    G4String PV_Name_metal;
    
    G4String Solid_Name_substrate;
    G4String LV_Name_substrate;
    G4String PV_Name_substrate;
    
    G4VisAttributes* TheVisAtt          = NULL;
    G4VisAttributes* TheVisAttMetal     = NULL;
    G4VisAttributes* TheVisAttSubstrate = NULL;
    
    if(ReadOut[sensorID] <= 0) {
      TString NameTmp;
      
      NameTmp    = TString("InnertMaterial_Solid") + long(sensorID+1);
      Solid_Name = NameTmp.Data();
      NameTmp    = TString("InnertMaterial_LV") + long(sensorID+1);
      LV_Name    = NameTmp.Data();
      NameTmp    = TString("InnertMaterial_PV") + long(sensorID+1);
      PV_Name    = NameTmp.Data();
      
      Some_LV = LV_Name;
      _ListOfAllLogicalVolume_Names.push_back(Some_LV);
      
      TheVisAtt = InnertVisAtt; 
    }
    else {
      TString NameTmp;
      
      //Sensitive layer
      NameTmp    = TString("MimosaSensor_Solid") + long(sensorID+1);
      Solid_Name = NameTmp.Data();
      NameTmp    = TString("MimosaSensor_LV") + long(sensorID+1);
      LV_Name    = NameTmp.Data();
      NameTmp    = TString("MimosaSensor_PV") + long(sensorID+1);
      PV_Name    = NameTmp.Data();
      
      _ListOfSensitiveDetectors_ID.push_back(sensorID+1);
      _ListOfSensitiveDetectors_Names.push_back(LV_Name);
      
      TheVisAtt = MimosaVisAtt;
      
      if(WithEpiMetalAndSubstrate) {
        //Metal layer
        NameTmp          = TString("MimosaSensorMetal_Solid") + long(sensorID+1);
        Solid_Name_metal = NameTmp.Data();
        NameTmp          = TString("MimosaSensorMetal_LV") + long(sensorID+1);
        LV_Name_metal    = NameTmp.Data();
        NameTmp          = TString("MimosaSensorMetal_PV") + long(sensorID+1);
        PV_Name_metal    = NameTmp.Data();
        _ListOfSensitiveDetectors_Names.push_back(LV_Name_metal);
        TheVisAttMetal     = MimosaMetalVisAtt;

        //Substrate layer
        NameTmp              = TString("MimosaSensorSubstrate_Solid") + long(sensorID+1);
        Solid_Name_substrate = NameTmp.Data();
        NameTmp              = TString("MimosaSensorSubstrate_LV") + long(sensorID+1);
        LV_Name_substrate    = NameTmp.Data();
        NameTmp              = TString("MimosaSensorSubstrate_PV") + long(sensorID+1);
        PV_Name_substrate    = NameTmp.Data();
        _ListOfSensitiveDetectors_Names.push_back(LV_Name_substrate);
        TheVisAttSubstrate   = MimosaSubstrateVisAtt;
      }
    }
    
    CheckMaterial(fMaterial[sensorID]);
    
    G4double thickness_fraction = 1.0/4.0;
    G4double MetalThickness     = 0.0;
    G4double EpiThickness       = fSensorThickness[sensorID];
    G4double SubstrateThickness = 0.0;
    if(WithEpiMetalAndSubstrate) {
      MetalThickness     = fSensorThickness[sensorID]*fSensorMetalThickness[sensorID];
      EpiThickness       = fSensorThickness[sensorID]*fSensorEpiThickness[sensorID];
      SubstrateThickness = fSensorThickness[sensorID]*(1.0 - fSensorEpiThickness[sensorID] - fSensorMetalThickness[sensorID]);
    }
    
    G4double maxStep = EpiThickness;
    if(WithEpiMetalAndSubstrate) {
      if(maxStep > MetalThickness)     maxStep = MetalThickness;
      if(maxStep > SubstrateThickness) maxStep = SubstrateThickness;
    }
    maxStep *= thickness_fraction;
    //fStepLimit = new G4UserLimits(maxStep);
    
    G4double MinimumThickness = 1.0e-3*um;
    
    if(EpiThickness < MinimumThickness) EpiThickness = MinimumThickness;
    
    //The sensitive voume
    G4Box* solidMimosa            = new G4Box(Solid_Name,fSensorLengthX[sensorID]/2,fSensorLengthY[sensorID]/2,EpiThickness/2);
    fLogicMimosaSensor[sensorID]  = new G4LogicalVolume(solidMimosa,_material_map[fMaterial[sensorID]],LV_Name);
    fLogicMimosaSensor[sensorID]->SetVisAttributes(TheVisAtt);
    if(ReadOut[sensorID] <= 0) _ListOfAllLogicalVolumes.push_back(fLogicMimosaSensor[sensorID]);
    fPhysicMimosaSensor[sensorID] = new G4PVPlacement(rot[sensorID],SensorPosition[sensorID],PV_Name,fLogicMimosaSensor[sensorID],worldPV, false,0,fCheckOverlaps);
    fLogicMimosaSensor[sensorID]->SetUserLimits(new G4UserLimits(maxStep));
    
    if(WithEpiMetalAndSubstrate) {
      if(MetalThickness     < MinimumThickness) MetalThickness     = MinimumThickness;
      if(SubstrateThickness < MinimumThickness) SubstrateThickness = MinimumThickness;
      G4ThreeVector delta;
      
      //The Metalization layer
      G4Box* solidMimosaMetal            = new G4Box(Solid_Name,fSensorLengthX[sensorID]/2,fSensorLengthY[sensorID]/2,MetalThickness/2);
      fLogicMimosaSensorMetal[sensorID]  = new G4LogicalVolume(solidMimosaMetal,_material_map[fMaterial[sensorID]],LV_Name_metal);
      fLogicMimosaSensorMetal[sensorID]->SetVisAttributes(TheVisAttMetal);
      _ListOfAllLogicalVolumes.push_back(fLogicMimosaSensorMetal[sensorID]);
      delta = G4ThreeVector(0.0,0.0,1.0);
      delta.transform(*Invrot[sensorID]);
      delta *= 0.5*(EpiThickness + MetalThickness)*(1.0 + 1.0e-4);
      G4ThreeVector SensorMetalPosition = SensorPosition[sensorID] + delta;
      fPhysicMimosaSensorMetal[sensorID] = new G4PVPlacement(rot[sensorID],SensorMetalPosition,PV_Name_metal,fLogicMimosaSensorMetal[sensorID],worldPV, false,0,fCheckOverlaps);
      fLogicMimosaSensorMetal[sensorID]->SetUserLimits(new G4UserLimits(maxStep));
      
      //The Substrate layer
      G4Box* solidMimosaSubstrate            = new G4Box(Solid_Name,fSensorLengthX[sensorID]/2,fSensorLengthY[sensorID]/2,SubstrateThickness/2);
      fLogicMimosaSensorSubstrate[sensorID]  = new G4LogicalVolume(solidMimosaSubstrate,_material_map[fMaterial[sensorID]],LV_Name_substrate);
      fLogicMimosaSensorSubstrate[sensorID]->SetVisAttributes(TheVisAttSubstrate);
      _ListOfAllLogicalVolumes.push_back(fLogicMimosaSensorSubstrate[sensorID]);
      delta = G4ThreeVector(0.0,0.0,1.0);
      delta.transform(*Invrot[sensorID]);
      delta *= 0.5*(EpiThickness + SubstrateThickness)*(1.0 + 1.0e-4);
      G4ThreeVector SensorSubstratePosition = SensorPosition[sensorID] - delta;
      fPhysicMimosaSensorSubstrate[sensorID] = new G4PVPlacement(rot[sensorID],SensorSubstratePosition,PV_Name_substrate,fLogicMimosaSensorSubstrate[sensorID],worldPV, false,0,fCheckOverlaps);
      fLogicMimosaSensorSubstrate[sensorID]->SetUserLimits(new G4UserLimits(maxStep));
    }
  }

  if(fSetup->CheckIfSourceExists(fSetup->TrackerParameter.ExpSetup)) {
    G4double maxStep = SourceHight/4.0;

    G4Tubs* solidSource  = new G4Tubs(Source_Solid,0.0,SourceRadius,SourceHight/2,0.0,2.0*TMath::Pi());
    fLogicSource         = new G4LogicalVolume(solidSource,_material_map[SourceMaterial],Source_LV);
    fLogicSource->SetVisAttributes(TheVisAttSource);
    _ListOfAllLogicalVolumes.push_back(fLogicSource);
    fPhysicSource        = new G4PVPlacement(Source_rot,SourcePosition,Source_PV,fLogicSource,worldPV, false,0,fCheckOverlaps);
    fLogicSource->SetUserLimits(new G4UserLimits(maxStep));
  }

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter
 
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void MimosaSimuDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  for(G4int iSensor = 0; iSensor < G4int(_ListOfSensitiveDetectors_ID.size()); iSensor++) {
    G4String MimosaSensorSDname = _ListOfSensitiveDetectors_Names[iSensor];
    G4cout << "Mimosa Sensor Sensitive Detector Name = " << MimosaSensorSDname.data() << G4endl;
    TString NameTmp = TString("MimosaSensorHitsCollection") + long(_ListOfSensitiveDetectors_ID[iSensor]);
    G4String MimosaSensorCollname = NameTmp.Data();
    cout << "MimosaSensitive collection name " << MimosaSensorCollname.data() << endl;
    
    MimosaSimuTrackerSD* aTrackerSD = new MimosaSimuTrackerSD(MimosaSensorSDname,
                                                              MimosaSensorCollname,
						              _ListOfSensitiveDetectors_ID[iSensor],
							      fSetup->GetAnalysisPar().SavePlots);
     SDman->AddNewDetector(aTrackerSD);
     fLogicMimosaSensor[_ListOfSensitiveDetectors_ID[iSensor]-1]->SetSensitiveDetector(aTrackerSD);
   }

   for(G4int iLV = 0; iLV < G4int(_ListOfAllLogicalVolume_Names.size()); iLV++) {
     G4String SDname = _ListOfAllLogicalVolume_Names[iLV];
     G4cout << "Trajectory Sensor Sensitive Detector Name = " << SDname.data() << G4endl;
     TString NameTmp = TString("TrajectoryHitsCollection") + long(iLV+1);
     G4String TrajectoryCollName = NameTmp.Data();
     cout << "Tracjectory collection name " << TrajectoryCollName.data() << endl;
     
     MimosaSimuTrajectorySD* aTrajectorySD = new MimosaSimuTrajectorySD(SDname,
									TrajectoryCollName,
							                fSetup->GetAnalysisPar().SavePlots);
     SDman->AddNewDetector(aTrajectorySD);
     _ListOfAllLogicalVolumes[iLV]->SetSensitiveDetector(aTrajectorySD);
   }

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = fSetup->TrackerParameter.BFieldDirection;
  fieldValue *= fSetup->TrackerParameter.BFieldMagnitude*tesla;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuDetectorConstruction::SetMimosaSensorMaterial(G4int sensorID, G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if(fSensorMaterial[sensorID] != pttoMaterial) {
    if( pttoMaterial ) {
      fSensorMaterial[sensorID] = pttoMaterial;
      if (fLogicMimosaSensor[sensorID]) fLogicMimosaSensor[sensorID]->SetMaterial(fChamberMaterial);
      G4cout << G4endl << "----> The sensor " << sensorID << " is made of " << materialName << G4endl;
    } else {
      G4cout << G4endl << "-->  WARNING from SetMimosaSensorMaterial : " << materialName << " not found" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MimosaSimuDetectorConstruction::CheckMaterial(TString Material)
{

  bool IsInList = false;
  for(int i=0;i<int(_material_name_list.size());i++) {
    if(Material == _material_name_list[i]) {
      IsInList = true;
      break;
    }
  }

  if(!IsInList) {
    G4cout << G4endl;
    G4cout << "MimosaSimuDetectorConstruction::The specified material " << Material.Data() << " is not in material list," << G4endl;
    for(int i=0;i<int(_material_name_list.size());i++) {
      G4cout << " - " << _material_name_list[i].Data() << G4endl;
    }
    G4cout << "MimosaSimuDetectorConstruction::Check your inputs. Exiting now!!!" << G4endl;
    G4cout << G4endl;
    
    assert(false);
  }
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  MimosaSimuDetectorConstruction::SetRotationMatrix(G4double RotationX,
							G4double RotationY,
							G4double RotationZ,
							G4RotationMatrix* &rot)
{
  
  if(fSetup->GetTrackerPar().DPrecAlignMethod == 0) {
    rot->rotateX(RotationX);
    rot->rotateY(RotationY);
    rot->rotateZ(RotationZ);
  }
  else if(fSetup->GetTrackerPar().DPrecAlignMethod == 1) {
    rot->rotateX(RotationX);
    rot->rotateY(RotationY);
    rot->rotateZ(RotationZ);
  }
  
  return;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......