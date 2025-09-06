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

    void AddScintillation_ev_depth(G4float z) {
      if ((z>-5) and (z<-4)) {ScintPhotons_0 += 1;}
      if ((z>-4) and (z<-3)) {ScintPhotons_1 += 1;}
      if ((z>-3) and (z<-2)) {ScintPhotons_2 += 1;}
      if ((z>-2) and (z<-1)) {ScintPhotons_3 += 1;}
      if ((z>-1) and (z<0)) {ScintPhotons_4 += 1;}
      if ((z>0) and (z<1)) {ScintPhotons_5 += 1;}
      if ((z>1) and (z<2)) {ScintPhotons_6 += 1;}
      if ((z>2) and (z<3)) {ScintPhotons_7 += 1;}
      if ((z>3) and (z<4)) {ScintPhotons_8 += 1;}
      if ((z>4) and (z<5)) {ScintPhotons_9 += 1;}

    }

// count mean scintillation depth
   void MeanScintDepth_ev(G4double z2) {
     Scint_depth = (Scint_depth * (ScintPhotons-1)+z2)/ScintPhotons;
     //Scint_depth_std = (Scint_depth_std * (ScintPhotons-1)+pow((Scint_depth-z2),2))/ScintPhotons;

    }


	G4double Absenergy;
	G4long AbsPhotonsSiPM;
	G4long AbsPhotonsScint;

	G4long AbsPhotonsSiPM_prior;
	G4long AbsPhotonsScint_prior;

	G4long ScintPhotons;
  G4long ScintPhotons_0;
  G4long ScintPhotons_1;
  G4long ScintPhotons_2;
  G4long ScintPhotons_3;
  G4long ScintPhotons_4;
  G4long ScintPhotons_5;
  G4long ScintPhotons_6;
  G4long ScintPhotons_7;
  G4long ScintPhotons_8;
  G4long ScintPhotons_9;
	G4double Scint_depth;
	G4double Scint_depth_std;

	G4long TotalInternalReflection;
	G4long FresnelRefraction;
	G4long FresnelReflection;
	G4long BoundAbsorption;

	G4long LambReflection;
	G4long SpikeReflection;

  private:
    RunAction*       runAction;
};
