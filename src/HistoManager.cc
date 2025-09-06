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
/// \file optical/OpNovice2/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id: HistoManager.cc 104417 2017-05-30 08:30:48Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("results")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetActivation(true);    // enable inactivation of histograms




//Ntuple for photons count.

  analysisManager->CreateNtuple("phot_count", "number of scintillated photons");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("counts");
  analysisManager->CreateNtupleIColumn("counts_0");
  analysisManager->CreateNtupleIColumn("counts_1");
  analysisManager->CreateNtupleIColumn("counts_2");
  analysisManager->CreateNtupleIColumn("counts_3");
  analysisManager->CreateNtupleIColumn("counts_4");
  analysisManager->CreateNtupleIColumn("counts_5");
  analysisManager->CreateNtupleIColumn("counts_6");
  analysisManager->CreateNtupleIColumn("counts_7");
  analysisManager->CreateNtupleIColumn("counts_8");
  analysisManager->CreateNtupleIColumn("counts_9");

  analysisManager->FinishNtuple();

  //analysisManager->CreateH2("absXY","XY distribution of photons absorbed in SiPM",
  //                          51, -25.5, 25.5, 51, -25.5, 25.5);
  analysisManager->CreateH2("absXY","XY distribution of photons absorbed in SiPM",
                            18, -27, 27, 17, -25.5, 25.5);

  analysisManager->CreateH2("X_ev","X distribution in event", 51, -25.5, 25.5, 1, 0, 1);
  analysisManager->CreateH2("R_ev","R distribution in event", 361, 0.1, 36.1, 100, 1, 100);

  analysisManager->CreateH1("sc_spec", "scintillation spectrum", 250, 1.77, 3.31);
  analysisManager->CreateH1("absR","R distribution of photons absorbed in SiPM",
                            3610, 0, 36.1);
  analysisManager->CreateH1("absX","X distribution of photons absorbed in SiPM",
                            510, -25.5, 25.5);

//Ntuple for absorption position.
  analysisManager->CreateNtuple("absorption", "absorption position XY");
  analysisManager->CreateNtupleIColumn(1, "eventID");
  analysisManager->CreateNtupleFColumn(1, "x");
  analysisManager->CreateNtupleFColumn(1, "y");
  analysisManager->FinishNtuple(1);

//Ntuple for scintillation spectrum.
  analysisManager->CreateNtuple("scintillation", "Scintillation spectrum");
	analysisManager->CreateNtupleFColumn(2, "energy");
	analysisManager->FinishNtuple(2);

//Ntuple for scintillation depth.
  analysisManager->CreateNtuple("scint_depth", "scintillation depth");
	analysisManager->CreateNtupleIColumn(3, "eventID");
	analysisManager->CreateNtupleFColumn(3, "z");
	analysisManager->FinishNtuple(3);

//Ntuple for absorption spectrum
	analysisManager->CreateNtuple("abs_sp", "absorption spectrum");
  analysisManager->CreateNtupleFColumn(4, "energy");
  analysisManager->FinishNtuple(4);

//Ntuple for opt. photons interaction
  analysisManager->CreateNtuple("status", "opt. photons interaction");
  analysisManager->CreateNtupleIColumn(5, "abs_SiPM");
  analysisManager->CreateNtupleIColumn(5, "abs_Scint");
  analysisManager->CreateNtupleIColumn(5, "abs_SiPM_prior");
  analysisManager->CreateNtupleIColumn(5, "abs_Scint_prior");
  analysisManager->CreateNtupleIColumn(5, "tot_int_ref");
  analysisManager->CreateNtupleIColumn(5, "fres_refrac");
  analysisManager->CreateNtupleIColumn(5, "fres_refl");
  analysisManager->CreateNtupleIColumn(5, "bound_abs");
  analysisManager->CreateNtupleIColumn(5, "lamb_refl");
  analysisManager->CreateNtupleIColumn(5, "spike_refl");
  analysisManager->FinishNtuple(5);

//Ntuple for primary interaction depth
  analysisManager->CreateNtuple("pr_int", "primary interaction depth");
	analysisManager->CreateNtupleIColumn(6, "eventID");
	analysisManager->CreateNtupleIColumn(6, "TrackID");
  analysisManager->CreateNtupleFColumn(6, "pr_int_depth");
	analysisManager->CreateNtupleIColumn(6, "interactionID");
  analysisManager->CreateNtupleFColumn(6, "energy_track");

  analysisManager->FinishNtuple(6);

  //Ntuple for scintillation depth.
  analysisManager->CreateNtuple("scint_depth_std", "scintillation depth std");
  analysisManager->CreateNtupleFColumn(7, "z");
  analysisManager->FinishNtuple(7);

  //Ntuple for phi and thetha.
  analysisManager->CreateNtuple("scat_angles", "scattering angle");
  analysisManager->CreateNtupleIColumn(8, "eventID");
  analysisManager->CreateNtupleDColumn(8, "phi_0");
  analysisManager->CreateNtupleDColumn(8, "phi");
  analysisManager->CreateNtupleDColumn(8, "theta_0");
  analysisManager->CreateNtupleDColumn(8, "theta");
  analysisManager->CreateNtupleDColumn(8, "Ekin_0");
  analysisManager->CreateNtupleDColumn(8, "Ekin_1");
  analysisManager->FinishNtuple(8);

}
