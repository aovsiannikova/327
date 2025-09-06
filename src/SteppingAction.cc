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
// $Id: SteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file optical/OpNovice2/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "TrackInformation.hh"
#include "Run.hh"

#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4ProcessManager.hh"

#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction* EvAct, G4bool primary_int)
:event(EvAct),
fPrimaryInt(primary_int)

{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4bool absorb_position=1; 	   //write down eventID, x (and y) for optical photons absorbed in SiPM 	nt_absorption.csv/h2_absXY.csv 	#1
  G4bool status         =0;		   //write down optical photons interaction		                   nt_status.csv		                #5
  G4bool abs_spectrum   =0;	     //spectrum of the optical photons absorbed in SiPM 		              nt_abs_sp.csv		                #4
  G4bool scint_spectrum =0;	     //spectrum of scintillated optical photons			                        nt_scintillation.csv	          #2
  G4bool scint_depth    =1;		   //z of the scintillation process			                    nt_scint_depth.csv	            #3
  theta         =false;		   //primary interaction of gamma photons                                 nt_pr_int.csv                   #6




  static G4ParticleDefinition* opticalphoton =
              G4OpticalPhoton::OpticalPhotonDefinition();

  static G4AnalysisManager* analysisMan = G4AnalysisManager::Instance();

  Run* run = static_cast<Run*>(
               G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4Track* track = step->GetTrack();
  G4int eventID 	= G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4StepPoint* endPoint   = step->GetPostStepPoint();
  G4StepPoint* startPoint = step->GetPreStepPoint();

  G4String particleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();

  TrackInformation* trackInfo =
                        (TrackInformation*)(track->GetUserInformation());

  if (particleName == "opticalphoton") {

    const G4VProcess* pds = endPoint->GetProcessDefinedStep();


    if (pds->GetProcessName() == "OpAbsorption") {
        run->AddOpAbsorption ();
        const G4String currentPhysicalName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

            if ((currentPhysicalName=="SiPMXY") or (currentPhysicalName=="SiPMYZ2") or (currentPhysicalName=="SiPMYZ1")
            or (currentPhysicalName=="SiPM_Box1") or (currentPhysicalName=="SiPM_Box2") or (currentPhysicalName=="SiPMTrap1")
            or (currentPhysicalName=="SiPMTrap2")) { //

                run->AddSiPMOpAbsorption();
                event->AddSiPMOpAbsorption_ev();
// *****************   get position of the optical photons absorbed in SiPM                ************************************************
                if (absorb_position) {
        			       G4ThreeVector dir = endPoint->GetPosition();
        			       G4double x1 = dir.x();
        			       G4double y1 = dir.y();
        			       // fill histogram with absorption position (XY)            ****************************************************
                     //analysisMan->FillH2(0, x1/mm, y1/mm);
                     // fill histogram with absorption position (R)            ****************************************************
                     //analysisMan->FillH1(1, sqrt(pow(x1,2) + pow(y1,2))/mm);
                     // fill histogram with absorption position (X)            ****************************************************
                     //analysisMan->FillH1(2, x1/mm);

                     // fill histogram with absorption position (X) and eventID            ****************************************************
                     //analysisMan->FillH2(1, x1/mm, eventID);
                      // fill histogram with absorption position (R) and eventID            ****************************************************
                     analysisMan->FillH2(2, sqrt(pow(x1,2) + pow(y1,2))/mm, eventID);
				             // fill ntuple with absorption position (XY)
				             //analysisMan->FillNtupleIColumn(1, 0, eventID);
                		 //analysisMan->FillNtupleFColumn(1, 1, x1/mm);
                		 //analysisMan->FillNtupleFColumn(1, 2, y1/mm);
        			       //analysisMan->AddNtupleRow(1);
                }
                if (abs_spectrum){

        			       G4double ekin = endPoint->GetKineticEnergy();
                                        analysisMan->FillNtupleFColumn(4, 0, ekin/eV);
        			       analysisMan->AddNtupleRow(4);
                }
            }
            if (currentPhysicalName=="Scintillator") {
                run->AddScOpAbsorption();
                event->AddScintOpAbsorption_ev();
            }

            if (trackInfo->GetIsFirstTankX()) {
                run->AddOpAbsorptionPrior();
                if (currentPhysicalName=="Scintillator") {
                    run->AddOpAbsorptionPriorScint();
                    event->AddOpAbsorptionPriorScint_ev();
                }
                if (currentPhysicalName=="SiPM") {
                    run->AddOpAbsorptionPriorSiPM();
                    event->AddOpAbsorptionPriorSiPM_ev();
                }
            }

    }
    else if (pds->GetProcessName() == "OpRayleigh") {
      run->AddRayleigh();
    }

    // optical process has endpt on bdry,
    if (status){
    	if (endPoint->GetStepStatus() == fGeomBoundary) {
      		const G4String NextPhysicalName = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      		const G4String currentPhysicalName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      		const G4DynamicParticle* theParticle = track->GetDynamicParticle();

          G4ThreeVector oldMomentumDir = theParticle->GetMomentumDirection();

	        G4ThreeVector m0 = startPoint->GetMomentumDirection();
      		G4ThreeVector m1 = endPoint->GetMomentumDirection();

          G4OpBoundaryProcessStatus theStatus = Undefined;

      		G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
      		G4int MAXofPostStepLoops =   OpManager->GetPostStepProcessVector()->entries();
      		G4ProcessVector* postStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

      		if (trackInfo->GetIsFirstTankX()) {

            if (NextPhysicalName!="SiPM") {
              trackInfo->SetIsFirstTankX(false);
            }
	       		run->AddTotalSurface();

	        	for (G4int i=0; i<MAXofPostStepLoops; ++i) {
		        	G4VProcess* currentProcess = (*postStepDoItVector)[i];

		        	G4OpBoundaryProcess* opProc = dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
          		if (opProc) {
                theStatus = opProc->GetStatus();

                if (theStatus == Transmission) {
                  G4cout << theStatus << G4endl;
            			run->AddTransmission();
                }

          			else if (theStatus == FresnelRefraction) {
                            if (NextPhysicalName!="SiPM") {
                                event->AddFresnelRefraction_ev();
                                run->AddFresnelRefraction();
                            }
                        }
			    		  else if (theStatus == FresnelReflection) {
                            event->AddFresnelReflection_ev();
                            run->AddFresnelReflection();
			    		  }
			    		  else if (theStatus == TotalInternalReflection) {
                            event->AddTotalInternalReflection_ev();
		      				run->AddTotalInternalReflection();
				        }
				        else if (theStatus == LambertianReflection) {
				            run->AddLambertianReflection();
				            event->AddLambertianReflection_ev();
                        }
				        else if (theStatus == LobeReflection) {
				            run->AddLobeReflection();
				        }
				        else if (theStatus == SpikeReflection) {
				            run->AddSpikeReflection();
				            event->AddSpikeReflection_ev();
                        }
			          else if (theStatus == BackScattering) {
				            run->AddBackScattering();
				        }
                else if (theStatus == Absorption) {
                            event->AddAbsorption_ev();
                            run->AddAbsorption();
                }
                else if (theStatus == SameMaterial) {
                            run->AddSameMaterial();
                }
		            else if (theStatus == StepTooSmall) {
				            run->AddStepTooSmall();
		            }
                else if (theStatus == NoRINDEX) {
	                 run->AddNoRINDEX();
		            }
                else {
                  G4cout << theStatus << G4endl;
                }

              } //if OpProc
            } //for
          } // if 1
      } //if ggem b
    } //if bound_pr
  } //if opt phot

  else { // particle != opticalphoton

		if ((particleName == "gamma") and (fPrimaryInt)) {

		    GammaPrimaryInteraction(step, analysisMan, eventID, run);

		}

    const std::vector<const G4Track*>* secondaries =
                                step->GetSecondaryInCurrentStep();
    for (auto sec : *secondaries) {
      	if (sec->GetDynamicParticle()->GetParticleDefinition() == opticalphoton){
            if (sec->GetCreatorProcess()
                    ->GetProcessName().compare("Scintillation") == 0) {
                //		Counts scintillation photons in a run.
                run->AddScintillation();


                if (scint_depth){
		                //*******************+ Get creation depth  *****************
                    G4ThreeVector vec = endPoint->GetPosition();
                    G4double z2 = vec.z()/mm;
                  // 		counts scintillation photons in event

                    G4float z_scint = sec->GetPosition().z()/mm;
                    event->AddScintillation_ev_depth(z2);

                  //event->AddScintillation_ev_depth2(z_scint);
	        	        // compute mean scintillation position                          ****************************************************
                    //event->MeanScintDepth_ev(z2);
                }
                // 		counts scintillation photons in event
                event->AddScintillation_ev();



                if (scint_spectrum){
                   //		Get emission spectra.               *************************************
                   G4double en = sec->GetKineticEnergy();
                   analysisMan->FillH1(0, en/eV);

                   //run->AddScintillationEnergy(en);

                   analysisMan->FillNtupleFColumn(2, 0, en/eV);
                   analysisMan->AddNtupleRow(2);
                }
            }
        }
    }
  }

  return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::GammaPrimaryInteraction(const G4Step* step, G4AnalysisManager* aMan, G4int evID, Run* run){
  G4StepPoint* end = step->GetPostStepPoint();
  G4StepPoint* start = step->GetPreStepPoint();

  G4ThreeVector vec = end->GetPosition();
  G4double z2 	= vec.z();
  G4Track* track = step->GetTrack();
  G4int trackID = track->GetTrackID();

  G4String proc =  end->GetProcessDefinedStep()->GetProcessName();
  G4float en_track 	= track->GetKineticEnergy();

  aMan->FillNtupleIColumn(6, 0, evID);
  aMan->FillNtupleIColumn(6, 1, trackID);
  aMan->FillNtupleFColumn(6, 2, z2/mm);

  if (proc =="conv"){
    // track->SetTrackStatus(fKillTrackAndSecondaries);
    run->AddConv();
    aMan->FillNtupleIColumn(6, 3, 0);		//conv - 0
    aMan->AddNtupleRow(6);
  }

  else if (proc =="compt"){
    run->AddCompt();
    aMan->FillNtupleIColumn(6, 3, 1);		//compt - 1
    aMan->AddNtupleRow(6);
    //  track->SetTrackStatus(fKillTrackAndSecondaries);

    if ((trackID==1) and (theta)) {
      G4ThreeVector momentum_gamma_0 = start->GetMomentum();
      double theta_gamma_0 = momentum_gamma_0.theta();
      double phi_gamma_0 = momentum_gamma_0.phi();
      double gamma_ekin_0 = start->GetKineticEnergy()/keV;
      G4ThreeVector momentum_gamma_1 = end->GetMomentum();
      double theta_gamma_1 = momentum_gamma_1.theta();
      double phi_gamma_1 = momentum_gamma_1.phi();
      double gamma_ekin_1 = end->GetKineticEnergy()/keV;

      aMan->FillNtupleIColumn(8, 0, evID);
      aMan->FillNtupleDColumn(8, 1, phi_gamma_0);
      aMan->FillNtupleDColumn(8, 2, phi_gamma_1);
      aMan->FillNtupleDColumn(8, 3, theta_gamma_0);
      aMan->FillNtupleDColumn(8, 4, theta_gamma_1);
      aMan->FillNtupleDColumn(8, 5, gamma_ekin_0);
      aMan->FillNtupleDColumn(8, 6, gamma_ekin_1);
      aMan->AddNtupleRow(8);
    }
  }

  else if (proc =="phot"){
    //track->SetTrackStatus(fKillTrackAndSecondaries);
    run->AddPhot();
    aMan->FillNtupleIColumn(6, 3, 2);		//phot - 2
    aMan->AddNtupleRow(6);
  }
  aMan->FillNtupleFColumn(6, 4, en_track/keV);

}