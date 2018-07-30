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
/// \file field/field01/src/TCEMfieldSetup.cc
/// \brief Implementation of the TCEMfieldSetup class
//
//
// $Id: TCEMfieldSetup.cc 77115 2013-11-21 15:06:37Z gcosmo $
//
//   User Field setup class implementation.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TCEMfieldSetup.hxx"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "G4UniformMagField.hh"
#include <Riostream.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TCEMfieldSetup::TCEMfieldSetup(TCEMfield* field)
 : fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fMagneticField(field),
   fStepper(0),
   fStepperType(0),
   fMinStep(0.)
{
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCEMfieldSetup::Initialize()
{
  fEquation       = new G4Mag_UsualEqRhs(fMagneticField);
  fMinStep        = 1.0*mm; // minimal step of 1 mm is default
  fStepperType    = 4;      // ClassicalRK4 is default stepper

  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  CreateStepperAndChordFinder();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCEMfieldSetup::~TCEMfieldSetup()
{
  delete fChordFinder;
  delete fStepper;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCEMfieldSetup::CreateStepperAndChordFinder()
{
  // Update field

  SetStepper();
  Info("CreateStepperAndChordFinder()", "The minimal step is equal to %f mm\n ", fMinStep/mm);

  fFieldManager->SetDetectorField(fMagneticField);
  fChordFinder = new G4ChordFinder( fMagneticField, fMinStep, fStepper);
  fFieldManager->SetChordFinder(fChordFinder);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCEMfieldSetup::SetStepper()
{
// Set stepper according to the stepper type

  if (fStepper) delete fStepper;

  switch (fStepperType)
  {
    case 0:
      fStepper = new G4ExplicitEuler(fEquation);
      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;
      break;
    case 1:
      fStepper = new G4ImplicitEuler(fEquation);
      G4cout<<"G4ImplicitEuler is called"<<G4endl;
      break;
    case 2:
      fStepper = new G4SimpleRunge(fEquation);
      G4cout<<"G4SimpleRunge is called"<<G4endl;
      break;
    case 3:
      fStepper = new G4SimpleHeum(fEquation);
      G4cout<<"G4SimpleHeum is called"<<G4endl;
      break;
    case 4:
      fStepper = new G4ClassicalRK4(fEquation, 8);
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
      break;
    case 5:
      fStepper = new G4HelixExplicitEuler(fEquation);
      G4cout<<"G4HelixExplicitEuler is called"<<G4endl;
      break;
    case 6:
      fStepper = new G4HelixImplicitEuler(fEquation);
      G4cout<<"G4HelixImplicitEuler is called"<<G4endl;
      break;
    case 7:
      fStepper = new G4HelixSimpleRunge(fEquation);
      G4cout<<"G4HelixSimpleRunge is called"<<G4endl;
      break;
    case 8:
      fStepper = new G4CashKarpRKF45(fEquation);
      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
      break;
    case 9:
      fStepper = new G4RKG3_Stepper(fEquation);
      G4cout<<"G4RKG3_Stepper is called"<<G4endl;
      break;
    default: fStepper = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
