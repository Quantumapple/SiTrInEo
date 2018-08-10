#ifndef MimosaSimuFieldSetup_h
#define MimosaSimuFieldSetup_h 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class MimosaSimuFieldMessenger;

class MimosaSimuFieldSetup
{
    public: 
      MimosaSimuFieldSetup();          // A zero field
      virtual ~MimosaSimuFieldSetup();

      void SetStepperType( G4int i ) { fStepperType = i; }

      void SetStepper();

      void SetMinStep(G4double s) { fMinStep = s; }

      void SetFieldValue(G4ThreeVector fieldVector);
      void SetFieldValue(G4double      fieldValue);
      G4ThreeVector GetConstantFieldValue();

      void UpdateField();

      G4FieldManager* GetLocalFieldManager() { return fLocalFieldManager; }

    protected:

      // Find the global Field Manager

      G4FieldManager*           GetGlobalFieldManager();

      G4FieldManager*           fFieldManager;
      G4FieldManager*           fLocalFieldManager;
      G4ChordFinder*            fChordFinder;
      G4ChordFinder*            fLocalChordFinder;
      G4Mag_UsualEqRhs*         fEquation;
      G4Mag_UsualEqRhs*         fLocalEquation;
      G4MagneticField*          fMagneticField;
      G4MagneticField*          fLocalMagneticField;

      G4MagIntegratorStepper*   fStepper;
      G4MagIntegratorStepper*   fLocalStepper;
      G4int                     fStepperType;

      G4double                  fMinStep;

      MimosaSimuFieldMessenger* fFieldMessenger;
};
#endif

