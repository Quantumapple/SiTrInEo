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
// $Id: MimosaSimu.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file MimosaSimu.cc
/// \brief Main program of the B2a example

#include "MimosaSimuActionInitialization.hh"
#include "MimosaSimuDetectorConstruction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UImanager.hh"
#include "MimosaSimuHistoManager.hh"
#include "MimosaSimuSetup.hh"
#include "TString.h"

#include "Randomize.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv) {
    // Detect interactive mode (if no arguments) and define UI session
    //

    // cout << "argc = " << argc << endl;

    TString TheConfigFile;
    G4UIExecutive *ui = 0;
    if (argc == 3) {
        int my_argc = 1;
        ui          = new G4UIExecutive(my_argc, argv);

        TheConfigFile = TString(argv[1]);
    } else if (argc == 2) {
        TheConfigFile = TString(argv[1]);
    } else {
        cout << endl;
        cout << "Need always to specify a config file as input" << endl;
        cout << endl;
        assert(false);
    }

    cout << endl;
    cout << "Specifying config file " << TheConfigFile.Data() << endl;
    cout << endl;

    MimosaSimuSetup *fSetup = new MimosaSimuSetup(TheConfigFile.Data());
    fSetup->ReadConfiguration();
    G4int verbosity = 1;
    if (fSetup->GetAnalysisPar().SavePlots)
        verbosity = 1;
    else
        verbosity = 0;

    // set an HistoManager
    MimosaSimuHistoManager *histo = new MimosaSimuHistoManager(fSetup);

    // Choose the Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

    // Construct the default run manager
    //
#ifdef G4MULTITHREADED
    // G4MTRunManager* runManager = new G4MTRunManager;
    G4RunManager *runManager = new G4RunManager;
#else
    G4RunManager *runManager = new G4RunManager;
#endif

    // Set mandatory initialization classes
    //
    runManager->SetUserInitialization(new MimosaSimuDetectorConstruction(fSetup));

    TString ThePhysicsList;
    // ThePhysicsList = TString("FTFP_BERT");
    // ThePhysicsList = TString("FTFP_BERT_HP");
    ThePhysicsList = TString("FTFP_BERT_LIV");
    // ThePhysicsList = TString("FTFP_BERT_EMV");
    // ThePhysicsList = TString("FTFP_BERT_EMX");
    // ThePhysicsList = TString("FTFP_BERT_EMZ");
    // ThePhysicsList = TString("FTFP_BERT_PEN);
    G4PhysListFactory factory;
    G4VModularPhysicsList *physicsList = factory.GetReferencePhysList(ThePhysicsList.Data());
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    // Set user action classes
    runManager->SetUserInitialization(new MimosaSimuActionInitialization(fSetup, histo));

    // Initialize visualization
    //
    G4VisManager *visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    //
    if (!ui) {
        // batch mode
        G4String command;
        TString TheGDMLCommand;

        command = "/run/initialize";
        UImanager->ApplyCommand(command);

        TheGDMLCommand = TString("/tracking/verbose  ") + long(verbosity);
        command        = TheGDMLCommand.Data();
        // command = "/tracking/verbose 0";
        UImanager->ApplyCommand(command);

        if (verbosity == 1)
            TheGDMLCommand = TString("/hits/verbose  2");
        else
            TheGDMLCommand = TString("/hits/verbose  0");
        // command = "/hits/verbose 0";
        // command = TheGDMLCommand.Data();
        UImanager->ApplyCommand(command);

        command = "/globalField/verbose 1";
        UImanager->ApplyCommand(command);

        TheGDMLCommand = TString("/run/beamOn  ") + long(fSetup->GetAnalysisPar().MCEvents);
        command        = TheGDMLCommand.Data();
        UImanager->ApplyCommand(command);

        // G4String command = "/control/execute ";
        // G4String fileName = argv[1];
        // G4String fileName = "run1.mac";
        // UImanager->ApplyCommand(command+fileName);
    } else {
        // interactive mode
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        if (ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute gui.mac");
        }
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !
    //

    delete visManager;
    cout << "MimosaSimu:: deleted visManager!!!" << endl;
    delete runManager;
    cout << "MimosaSimu:: deleted runManager!!!" << endl;
    delete histo;
    cout << "MimosaSimu:: deleted histo!!!" << endl;
    delete fSetup;
    cout << "MimosaSimu:: deleted fSetup!!!" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
