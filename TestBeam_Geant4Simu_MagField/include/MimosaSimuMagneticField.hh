#ifndef __MimosaSimuMagneticField_h__
#define __MimosaSimuMagneticField_h__

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

class MimosaSimuMagneticField : public G4MagneticField {
  public:
    MimosaSimuMagneticField();
    ~MimosaSimuMagneticField(){};

    inline void SetDirection(const G4ThreeVector &a) { fBDir = a; };
    inline G4ThreeVector GetDirection() const { return fBDir; }

    inline void SetMagnitude(const G4double a) { fBMag = a; };
    inline G4double GetMagnitude() const { return fBMag; }
    virtual void GetFieldValue(const G4double[4], G4double *) const;

  private:
    G4ThreeVector fBDir;
    G4double fBMag;
};

#endif
