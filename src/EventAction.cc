


#include "EventAction.hh"
#include "G4Event.hh"
#include "RunAction.hh"
#include "G4UserEventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "globals.hh"

using namespace std;
using namespace CLHEP;

EventAction::EventAction(RunAction* )
:G4UserEventAction()
{

}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{

	Absenergy = 0.;
	AbsPhotonsSiPM = 0;
	AbsPhotonsScint = 0;

	AbsPhotonsSiPM_prior = 0;
	AbsPhotonsScint_prior = 0;

	TotalInternalReflection = 0;
	FresnelRefraction	= 0;

	FresnelReflection	= 0;
	BoundAbsorption		= 0;

	LambReflection		= 0;
	SpikeReflection		= 0;

	ScintPhotons = 0;
  ScintPhotons_0 = 0;
  ScintPhotons_1 = 0;
  ScintPhotons_2 = 0;
  ScintPhotons_3 = 0;
  ScintPhotons_4 = 0;
  ScintPhotons_5 = 0;
  ScintPhotons_6 = 0;
  ScintPhotons_7 = 0;
  ScintPhotons_8 = 0;
  ScintPhotons_9 = 0;

	Scint_depth = 0.;
  Scint_depth_std = 0.;

//	G4cout << "Begin of event" << G4endl;
}

void EventAction::EndOfEventAction(const G4Event* evt) {

  // Accumulate statistics
  //
  G4int eventID = evt->GetEventID();

  // get analysis manager. add "1" if there is another Ntuple

  if (ScintPhotons > 0) {
    static G4AnalysisManager *analysisMan = G4AnalysisManager::Instance();
    analysisMan->FillNtupleIColumn(0, eventID);
    analysisMan->FillNtupleIColumn(1, ScintPhotons);
    analysisMan->FillNtupleIColumn(2, ScintPhotons_0);
    analysisMan->FillNtupleIColumn(3, ScintPhotons_1);
    analysisMan->FillNtupleIColumn(4, ScintPhotons_2);
    analysisMan->FillNtupleIColumn(5, ScintPhotons_3);
    analysisMan->FillNtupleIColumn(6, ScintPhotons_4);
    analysisMan->FillNtupleIColumn(7, ScintPhotons_5);
    analysisMan->FillNtupleIColumn(8, ScintPhotons_6);
    analysisMan->FillNtupleIColumn(9, ScintPhotons_7);
    analysisMan->FillNtupleIColumn(10, ScintPhotons_8);
    analysisMan->FillNtupleIColumn(11, ScintPhotons_9);

    analysisMan->AddNtupleRow(0);

    analysisMan->FillNtupleIColumn(3, 0, eventID);
    analysisMan->FillNtupleFColumn(3, 1, Scint_depth / mm);
    analysisMan->AddNtupleRow(3);


    analysisMan->FillNtupleIColumn(5, 0, AbsPhotonsSiPM);
    analysisMan->FillNtupleIColumn(5, 1, AbsPhotonsScint);
    analysisMan->FillNtupleIColumn(5, 2, AbsPhotonsSiPM_prior);
    analysisMan->FillNtupleIColumn(5, 3, AbsPhotonsScint_prior);
    analysisMan->FillNtupleIColumn(5, 4, TotalInternalReflection);
    analysisMan->FillNtupleIColumn(5, 5, FresnelRefraction);
    analysisMan->FillNtupleIColumn(5, 6, FresnelReflection);
    analysisMan->FillNtupleIColumn(5, 7, BoundAbsorption);
    analysisMan->FillNtupleIColumn(5, 8, LambReflection);
    analysisMan->FillNtupleIColumn(5, 9, SpikeReflection);
    analysisMan->AddNtupleRow(5);
  }

//  analysisMan->FillNtupleFColumn(7, 0, pow(Scint_depth_std,0.5)/mm);
//  analysisMan->AddNtupleRow(7);

  // Print per event (modulo n)
  //
//  G4int eventID = 1 + evt->GetEventID();
 //G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
//  G4int printModulo = 100;
//  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) )
//  {
//    G4cout << "---> End of event: " << eventID << G4endl;
//G4cout << "End of event " << eventID << " .  " << AbsPhotons << " photons absorbed in SiPM, " << ScintPhotons << " optical photons produced in scintillator." <<G4endl;
//  }
}
