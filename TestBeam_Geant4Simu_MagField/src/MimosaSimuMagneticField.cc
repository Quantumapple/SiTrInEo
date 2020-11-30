#include "MimosaSimuMagneticField.hh"

MimosaSimuMagneticField::MimosaSimuMagneticField() : fBDir(), fBMag(0.){};

void MimosaSimuMagneticField::GetFieldValue(const G4double nowPos[4], G4double *resField) const {
    G4ThreeVector lNowMag = fBMag * fBDir;
    resField[0]           = lNowMag[0];
    resField[1]           = lNowMag[1];
    resField[2]           = lNowMag[2];
};
