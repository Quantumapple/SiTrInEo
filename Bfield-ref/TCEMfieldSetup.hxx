//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TCEMfieldSetup_h
#define TCEMfieldSetup_h 1

#include "G4MagneticField.hh"
#include "TCEMfield.hxx"

#include "TVector3.h"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class TCEMfieldSetupMessenger;
class G4MagneticField;

class TCEMfieldSetup
{
public:
   TCEMfieldSetup(TCEMfield* field);
  virtual ~TCEMfieldSetup();

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
  TCEMfield*              fMagneticField;

  G4MagIntegratorStepper* fStepper;
  Int_t                   fStepperType;

  Double_t                fMinStep;
 
};

#endif
