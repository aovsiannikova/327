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
/// \file optical/OpNovice2/src/Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <numeric>

#include "Run.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::Run()
: G4Run()
{
  fParticle = nullptr;
  fEkin = -1.;

  fConv =0;
  fCompt =0;
  fPhot =0;

  fOtherGammaInt=0;

  fScintEnergy = 0.0;

  fScintCount = 0;
  fRayleighCount = 0;

  fOpAbsorption = 0;
  fSiPMOpAbsorption = 0;
  fScOpAbsorption = 0;
  fOpAbsorptionPrior = 0;
  fOpAbsorptionPriorSiPM = 0;
  fOpAbsorptionPriorScint = 0;

  fTotalSurface = 0;

  fBoundaryProcs.clear();
  fBoundaryProcs.resize(40);
  for (G4int i = 0; i < 40; ++i) {
    fBoundaryProcs[i] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
  fParticle = particle;
  fEkin = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  fConv 	+=localRun->fConv;
  fCompt 	+=localRun->fCompt;
  fPhot 	+=localRun->fPhot;

  fOtherGammaInt  +=localRun->fOtherGammaInt;

  fScintEnergy    += localRun->fScintEnergy;

  fScintCount     += localRun->fScintCount;
  fRayleighCount  += localRun->fRayleighCount;

  fTotalSurface   += localRun->fTotalSurface;

  fOpAbsorption   += localRun->fOpAbsorption;

// Photons absorbed in SiPM
  fSiPMOpAbsorption   += localRun->fSiPMOpAbsorption;

// in Scintillator
  fScOpAbsorption += localRun->fScOpAbsorption;

  fOpAbsorptionPrior += localRun->fOpAbsorptionPrior;
// absorbed in SiPM before bndr proc
  fOpAbsorptionPriorSiPM += localRun->fOpAbsorptionPriorSiPM;
// absorbed in Scintillator before bndr proc
  fOpAbsorptionPriorScint += localRun->fOpAbsorptionPriorScint;


  for (size_t i = 0; i < fBoundaryProcs.size(); ++i) {
    fBoundaryProcs[i] += localRun->fBoundaryProcs[i];
  }

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Run::EndOfRun()
{
  G4int TotNbofEvents = numberOfEvent;
  if (TotNbofEvents == 0) return;

  std::ios::fmtflags mode = G4cout.flags();
  G4int prec = G4cout.precision(2);

  G4cout << "\n    Run Summary\n";
  G4cout <<   "---------------------------------\n";
  G4cout << "Primary particle was: " << fParticle->GetParticleName()
         << " with energy " << G4BestUnit(fEkin/MeV, "Energy") << "." << G4endl;

  G4cout << "Primary interactions in a run: " << G4endl;
  if (fConv>0){
	G4cout << "pair productions: " << fConv << G4endl;
  }
  if (fCompt>0){
	G4cout << "Compton scattering: " << fCompt << G4endl;
  }
  if (fPhot>0){
	G4cout << "Photoelectric effect: " << fPhot << G4endl;
  }

  if (fOtherGammaInt>0){
	G4cout << "others: " << fOtherGammaInt << G4endl;
  }

  if (fParticle->GetParticleName() != "opticalphoton") {
    G4cout << "Total number of scintillation photons created: "
           << fScintCount << G4endl;
    if (fScintCount > 0) {
      G4cout << " Average energy: " << (fScintEnergy/eV)/fScintCount << " eV."
             << G4endl;
    }
  }
  G4cout.precision(5);
  G4cout << "OpAbsorption:      " << fOpAbsorption << "                of total created, %: " << fOpAbsorption*100/fScintCount   << G4endl;

// Photons absorbed in SiPM
  G4cout << "OpAbsorption per event in SiPM:      " << fSiPMOpAbsorption         << G4endl;

// Photons absorbed in Scintillator
  G4cout << "OpAbsorption per event in Scintillator:      " << fScOpAbsorption   << G4endl;



  G4cout << "\nSurface events this run:" << G4endl;
  G4cout << "# of primary particles:      " << std::setw(8) << TotNbofEvents
         << G4endl;
  G4cout << "OpAbsorption before surface: " << std::setw(8)
         << fOpAbsorptionPrior << G4endl;
  G4cout << "OpAbsorption before surface in SiPM: " << std::setw(8)
         << fOpAbsorptionPriorSiPM << G4endl;
  G4cout << "OpAbsorption before surface in scintillator: " << std::setw(8)
         << fOpAbsorptionPriorScint << G4endl;

  G4cout << "Total # of surface events:   " << std::setw(8) << fTotalSurface
         << G4endl;
  if (fParticle->GetParticleName() == "opticalphoton") {
    G4cout << "Unaccounted for:             " << std::setw(8)
         << fTotalSurface + fOpAbsorptionPrior - TotNbofEvents << G4endl;
  }
  G4cout << "\nSurface events by process:" << G4endl;
  if (fBoundaryProcs[Transmission] > 0) {
    G4cout << "  Transmission:              " << std::setw(8)
           << fBoundaryProcs[Transmission] << G4endl;
  }
  if (fBoundaryProcs[FresnelRefraction] > 0) {
    G4cout << "  Fresnel refraction:        " << std::setw(8)
           << fBoundaryProcs[FresnelRefraction] << G4endl;
  }
  if (fBoundaryProcs[FresnelReflection] > 0) {
    G4cout << "  Fresnel reflection:        " << std::setw(8)
           << fBoundaryProcs[FresnelReflection] << G4endl;
  }
  if (fBoundaryProcs[TotalInternalReflection] > 0) {
    G4cout << "  Total internal reflection: " << std::setw(8)
           << fBoundaryProcs[TotalInternalReflection] << G4endl;
  }
  if (fBoundaryProcs[LambertianReflection] > 0) {
    G4cout << "  Lambertian reflection:     " << std::setw(8)
           << fBoundaryProcs[LambertianReflection] << G4endl;
  }
  if (fBoundaryProcs[LobeReflection] > 0) {
    G4cout << "  Lobe reflection:           " << std::setw(8)
           << fBoundaryProcs[LobeReflection] << G4endl;
  }
  if (fBoundaryProcs[SpikeReflection] > 0) {
    G4cout << "  Spike reflection:          " << std::setw(8)
           << fBoundaryProcs[SpikeReflection] << G4endl;
  }
  if (fBoundaryProcs[BackScattering] > 0) {
    G4cout << "  Backscattering:            " << std::setw(8)
           << fBoundaryProcs[BackScattering] << G4endl;
  }
  if (fBoundaryProcs[Absorption] > 0) {
    G4cout << "  Absorption:                " << std::setw(8)
           << fBoundaryProcs[Absorption] << G4endl;
  }
  if (fBoundaryProcs[Detection] > 0) {
    G4cout << "  Detection:                 " << std::setw(8)
           << fBoundaryProcs[Detection] << G4endl;
  }
  if (fBoundaryProcs[NotAtBoundary] > 0) {
    G4cout << "  Not at boundary:           " << std::setw(8)
           << fBoundaryProcs[NotAtBoundary] << G4endl;
  }
  if (fBoundaryProcs[SameMaterial] > 0) {
    G4cout << "  Same material:             " << std::setw(8)
           << fBoundaryProcs[SameMaterial] << G4endl;
  }
  if (fBoundaryProcs[StepTooSmall] > 0) {
    G4cout << "  Step too small:            " << std::setw(8)
           << fBoundaryProcs[StepTooSmall] << G4endl;
  }
  if (fBoundaryProcs[NoRINDEX] > 0) {
    G4cout << "  No RINDEX:                 " << std::setw(8)
           << fBoundaryProcs[NoRINDEX] << G4endl;
  }
  // LBNL polished
  if (fBoundaryProcs[PolishedLumirrorAirReflection] > 0) {
    G4cout << "  Polished Lumirror Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedLumirrorAirReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedLumirrorGlueReflection] > 0) {
    G4cout << "  Polished Lumirror Glue reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedLumirrorGlueReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedAirReflection] > 0) {
    G4cout << "  Polished Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedAirReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedTeflonAirReflection] > 0) {
    G4cout << "  Polished Teflon Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedTeflonAirReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedTiOAirReflection] > 0) {
    G4cout << "  Polished TiO Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedTiOAirReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedTyvekAirReflection] > 0) {
    G4cout << "  Polished Tyvek Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedTyvekAirReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedVM2000AirReflection] > 0) {
    G4cout << "  Polished VM2000 Air reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedVM2000AirReflection] << G4endl;
  }
  if (fBoundaryProcs[PolishedVM2000GlueReflection] > 0) {
    G4cout << "  Polished VM2000 Glue reflection: " << std::setw(8)
           << fBoundaryProcs[PolishedVM2000GlueReflection] << G4endl;
  }
  // LBNL etched
  if (fBoundaryProcs[EtchedLumirrorAirReflection] > 0) {
    G4cout << "  Etched Lumirror Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedLumirrorAirReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedLumirrorGlueReflection] > 0) {
    G4cout << "  Etched Lumirror Glue reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedLumirrorGlueReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedAirReflection] > 0) {
    G4cout << "  Etched Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedAirReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedTeflonAirReflection] > 0) {
    G4cout << "  Etched Teflon Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedTeflonAirReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedTiOAirReflection] > 0) {
    G4cout << "  Etched TiO Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedTiOAirReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedTyvekAirReflection] > 0) {
    G4cout << "  Etched Tyvek Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedTyvekAirReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedVM2000AirReflection] > 0) {
    G4cout << "  Etched VM2000 Air reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedVM2000AirReflection] << G4endl;
  }
  if (fBoundaryProcs[EtchedVM2000GlueReflection] > 0) {
    G4cout << "  Etched VM2000 Glue reflection: " << std::setw(8)
           << fBoundaryProcs[EtchedVM2000GlueReflection] << G4endl;
  }
  // LBNL ground
  if (fBoundaryProcs[GroundLumirrorAirReflection] > 0) {
    G4cout << "  Ground Lumirror Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundLumirrorAirReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundLumirrorGlueReflection] > 0) {
    G4cout << "  Ground Lumirror Glue reflection: " << std::setw(8)
           << fBoundaryProcs[GroundLumirrorGlueReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundAirReflection] > 0) {
    G4cout << "  Ground Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundAirReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundTeflonAirReflection] > 0) {
    G4cout << "  Ground Teflon Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundTeflonAirReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundTiOAirReflection] > 0) {
    G4cout << "  Ground TiO Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundTiOAirReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundTyvekAirReflection] > 0) {
    G4cout << "  Ground Tyvek Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundTyvekAirReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundVM2000AirReflection] > 0) {
    G4cout << "  Ground VM2000 Air reflection: " << std::setw(8)
           << fBoundaryProcs[GroundVM2000AirReflection] << G4endl;
  }
  if (fBoundaryProcs[GroundVM2000GlueReflection] > 0) {
    G4cout << "  Ground VM2000 Glue reflection: " << std::setw(8)
           << fBoundaryProcs[GroundVM2000GlueReflection] << G4endl;
  }

  G4int sum = std::accumulate(fBoundaryProcs.begin(), fBoundaryProcs.end(), 0);
  G4cout << " Sum:                        " << std::setw(8) << sum << G4endl;
  G4cout << " Unaccounted for:            " << std::setw(8)
         << fTotalSurface - sum << G4endl;

  G4cout <<   "---------------------------------\n";

  G4cout.setf(mode, std::ios::floatfield);
  G4cout.precision(prec);
}
