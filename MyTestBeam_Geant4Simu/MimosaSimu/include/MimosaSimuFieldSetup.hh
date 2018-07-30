//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MimosaSimuFieldSetup_h
#define MimosaSimuFieldSetup_h 1

#include "G4MagneticField.hh"
#include "MimosaSimuField.hh"

#include "TVector3.h"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class MimosaSimuFieldSetupMessenger;
class G4MagneticField;

class MimosaSimuFieldSetup
{
public:
  MimosaSimuFieldSetup();
  MimosaSimuFieldSetup(MimosaSimuField* field);
  virtual ~MimosaSimuFieldSetup();

  void SetStepperType(Int_t i) { fStepperType = i; }
  void SetMinStep(Double_t s)  { fMinStep = s;     }
  
  void SetStepper();
  void Initialize();    //  Set parameters and call method below
  void CreateStepperAndChordFinder();
   
  G4FieldManager*  GetFieldManager() { return fFieldManager; }

protected:
  G4FieldManager*         fFieldManager;
  G4ChordFinder*          fChordFinder;
  G4Mag_UsualEqRhs*       fEquation;
  MimosaSimuField*        fMagneticField;

  G4MagIntegratorStepper* fStepper;
  Int_t                   fStepperType;

  Double_t                fMinStep;
 
};

#endif




