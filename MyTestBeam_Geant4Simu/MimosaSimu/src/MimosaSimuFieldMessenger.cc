#include "MimosaSimuFieldMessenger.hh"

#include "MimosaSimuFieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4SystemOfUnits.hh"

MimosaSimuFieldMessenger::MimosaSimuFieldMessenger(MimosaSimuFieldSetup* fieldSetup)
 : G4UImessenger(),
   fEMfieldSetup(fieldSetup),
   fFieldDir(0),
   fStepperCmd(0),
   fMagFieldCmd(0),
   fMinStepCmd(0),
   fUpdateCmd(0)
{
  fFieldDir = new G4UIdirectory("/field/");
  fFieldDir->SetGuidance("MimosaSimu field tracking control.");

  fStepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  fStepperCmd->SetGuidance("Select stepper type for magnetic field");
  fStepperCmd->SetParameterName("choice",true);
  fStepperCmd->SetDefaultValue(4);
  fStepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
 
  fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/field/setFieldZ",this);
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false,false);
  fMagFieldCmd->SetDefaultUnit("tesla");
  fMagFieldCmd->AvailableForStates(G4State_Idle);
 
  fMinStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);
  fMinStepCmd->SetGuidance("Define minimal step");
  fMinStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMinStepCmd->SetParameterName("min step",false,false);
  fMinStepCmd->SetDefaultUnit("mm");
  fMinStepCmd->AvailableForStates(G4State_Idle);
}

MimosaSimuFieldMessenger::~MimosaSimuFieldMessenger()
{
  delete fStepperCmd;
  delete fMagFieldCmd;
  delete fMinStepCmd;
  delete fFieldDir;
  delete fUpdateCmd;
}

void MimosaSimuFieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{
  if( command == fStepperCmd )
    fEMfieldSetup->SetStepperType(fStepperCmd->GetNewIntValue(newValue));
  if( command == fUpdateCmd )
    fEMfieldSetup->UpdateField();
  if( command == fMagFieldCmd )
  {
    fEMfieldSetup->SetFieldValue(fMagFieldCmd->GetNewDoubleValue(newValue));
    // Check the value
    G4cout << "Set field value to " <<
      fEMfieldSetup->GetConstantFieldValue() / gauss << " Gauss " << G4endl;
  }
  if( command == fMinStepCmd )
    fEMfieldSetup->SetMinStep(fMinStepCmd->GetNewDoubleValue(newValue));
}


