#ifndef MimosaSimuFieldMessenger_h
#define MimosaSimuFieldMessenger_h 1

#include "G4UImessenger.hh"

class MimosaSimuFieldSetup;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class MimosaSimuFieldMessenger: public G4UImessenger
{
    public:
      MimosaSimuFieldMessenger(MimosaSimuFieldSetup* );
      virtual ~MimosaSimuFieldMessenger();

      virtual void SetNewValue(G4UIcommand*, G4String);

    private:

      MimosaSimuFieldSetup*       fEMfieldSetup;

      G4UIdirectory*              fFieldDir;
      
      G4UIcmdWithAnInteger*       fStepperCmd;
      G4UIcmdWithADoubleAndUnit*  fMagFieldCmd;
      G4UIcmdWithADoubleAndUnit*  fMinStepCmd;
      G4UIcmdWithoutParameter*    fUpdateCmd;
};

#endif
