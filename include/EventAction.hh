#include "G4UserEventAction.hh"
#include "RunAction.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "HistoManager.hh"

class RunAction;
class G4Event;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction*);
    ~EventAction();


  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

 // count photons absorbed in SiPM
    void AddSiPMOpAbsorption_ev(void) {AbsPhotonsSiPM += 1;}
 // count photons absorbed in scintillator
    void AddScintOpAbsorption_ev(void) {AbsPhotonsScint += 1;}

 // count photons absorbed in SiPM directly
    void AddOpAbsorptionPriorSiPM_ev(void) {AbsPhotonsSiPM_prior += 1;}
 // count photons absorbed in scintillator directly
    void AddOpAbsorptionPriorScint_ev(void) {AbsPhotonsScint_prior += 1;}

 // count total internal reflection
    void AddTotalInternalReflection_ev(void) {TotalInternalReflection += 1;}

 // count Fresnel refractions (not on CeBr3 SiPM boundary)
    void AddFresnelRefraction_ev(void) {FresnelRefraction += 1;}

 // count Fresnel reflection
    void AddFresnelReflection_ev(void) {FresnelReflection += 1;}

 // count absorption on the boundary
    void AddAbsorption_ev(void) {BoundAbsorption += 1;}

 // count lambertian reflection
    void AddLambertianReflection_ev(void) {LambReflection += 1;}

 // count specular spike reflection
    void AddSpikeReflection_ev(void) {SpikeReflection += 1;}

 // count scintillation photons
    void AddScintillation_ev(void) {ScintPhotons += 1;}


	G4double Absenergy;
	G4int AbsPhotonsSiPM;
	G4int AbsPhotonsScint;

	G4int AbsPhotonsSiPM_prior;
	G4int AbsPhotonsScint_prior;

	G4int ScintPhotons;
	G4int TotalInternalReflection;
	G4int FresnelRefraction;
	G4int FresnelReflection;
	G4int BoundAbsorption;

	G4int LambReflection;
	G4int SpikeReflection;

  private:
    RunAction*       runAction;
};
