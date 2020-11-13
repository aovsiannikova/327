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
/// \file optical/OpNovice2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fDetectorMessenger(nullptr)
{
  fExpHall_x = fExpHall_y = fExpHall_z = 10.0*cm;
  fTank_x    = fTank_y    = fTank_z    =  5.0*cm;

  isQuadratic = 1;
  hasReflector=0;
  Reflector	=3;
  //0 - A, only diff reflection with 0.97 reflectivity;
  //1 - B, Spec Spike reflection on scint-air surface (sigalpha=0), diff reflection on air-tefl Surface;
  //2 - C, Spec Spike reflection on scint-air surface, Spec Spike on air-tefl surface;
  //3 - D, diffuse reflection on scint-air surface (sigalpha = 0.6), diffuse reflection on air-tefl surface.

  if (isQuadratic) {
    //fScintillator_a    =  35*mm; // 70 x 70 mm quadratic crystal a=70mm/2 = 35 mm
   fScintillator_a    =  25.5*mm;  // 51 x 51 mm quadratic crystal a=51mm/2 = 25.5 mm
  //fScintillator_a    =  40*mm;  //
  }

  else {
   fScintillator_r = 12.2 *mm; // 2.54cm = 25.4 mm diameter of cylindrical crystal r=25.4/2=12.2 mm
 //  fScintillator_r = 22.5 *mm; // 51 mm diameter of cylindrical crystal          r=51mm/2 = 25.5 mm - radius
  }
  fSiPM_z = 1.5*mm;

  //fScintillator_z       =  5*mm;  // 1 cm thick crystal
  fScintillator_z       =  10*mm;  // 2 cm thick crystal
  //fScintillator_z       =  12.7*mm;  // 2.54 cm thick crystal


  fTank = nullptr;
  fSurface = nullptr;


  fWorldMPT = new G4MaterialPropertiesTable();
  fScintillatorMPT = new G4MaterialPropertiesTable();
  fReflectorMPT = new G4MaterialPropertiesTable();
  fSiPMMPT = new G4MaterialPropertiesTable();

  fSurfaceMPT   = new G4MaterialPropertiesTable();
  fSc_RefMPT    = new G4MaterialPropertiesTable();
  fSc_SiPMMPT   = new G4MaterialPropertiesTable();
  fSi_RefMPT    = new G4MaterialPropertiesTable();


  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{



  const G4int nEntries = 25;


  G4double PhotonEnergy[nEntries] = {2.86*eV, 2.92*eV, 2.95*eV, 2.97*eV,
                                     3.00*eV, 3.03*eV, 3.06*eV, 3.09*eV, 3.12*eV,
                                     3.16*eV, 3.19*eV, 3.22*eV, 3.25*eV, 3.29*eV,
                                     3.32*eV, 3.35*eV, 3.36*eV, 3.40*eV, 3.42*eV,
                                     3.43*eV, 3.44*eV, 3.48*eV, 3.52*eV, 3.56*eV, 3.60*eV};


// ------------- Materials -------------

  G4NistManager* man = G4NistManager::Instance();


  G4Material* Si  =     man->FindOrBuildMaterial("G4_Si");
  G4Material* Al  =     man->FindOrBuildMaterial("G4_Al");


// The elements we need to build complexe materials

  G4Element* elCe   = man->FindOrBuildElement("Ce");
  G4Element* elBr   = man->FindOrBuildElement("Br");
  G4Element* elC    = man->FindOrBuildElement("C");
  G4Element* elF    = man->FindOrBuildElement("F");


// CeBr3 (the scintillator material)

    CeBr_3 = new G4Material(name = "CeBr_3", density = 5.2*g/cm3, ncomponents = 2);
    CeBr_3->AddElement(elCe, fractionmass = 0.36888);
    CeBr_3->AddElement(elBr, fractionmass = 0.63112);

// vacuum - tank, world
    Vacuum = new G4Material(name = "Vacuum", z = 1., a = 1.008*g/mole, density = 1.e-25*g/cm3);


// Teflon (the reflector material)

    C2F4 = new G4Material(name = "C2F4", density = 2.2*g/cm3, ncomponents = 2);
    C2F4->AddElement(elC, fractionmass = 0.24);
    C2F4->AddElement(elF, fractionmass = 0.76);



// ------------- Optical properties -------------

// CeBr3 optical properties

// relative scintillator yield
    G4double FastComp[nEntries] =     {0.010, 0.031, 0.061, 0.122,
                                    0.184, 0.245, 0.408, 0.551, 0.673,
                                    0.796, 0.816, 0.837, 0.857, 0.898,
                                    0.980, 1.000, 0.959, 0.816, 0.653,
                                    0.571, 0.408, 0.224, 0.143, 0.082, 0.020};
// Absorption length = Attenuation length
    G4double Absorption[nEntries] =   {155*mm, 150*mm, 147*mm, 145*mm,
                                     143*mm, 141*mm, 140*mm, 135*mm, 133*mm,
                                     130*mm, 129*mm, 126*mm, 125*mm, 121*mm,
                                     120*mm, 117*mm, 116*mm, 115*mm, 114*mm,
                                     113*mm, 112*mm, 110*mm, 108*mm, 106*mm, 105*mm};


    G4double rIndex[nEntries] =       {1.9, 1.9, 1.9, 1.9,
                                    1.9, 1.9, 1.9, 1.9, 1.9,
                                    1.9, 1.9, 1.9, 1.9, 1.9,
                                    1.9, 1.9, 1.9, 1.9, 1.9,
                                    1.9, 1.9, 1.9, 1.9, 1.9, 1.9};

    fScintillatorMPT->AddProperty("FASTCOMPONENT",PhotonEnergy, FastComp, nEntries)->SetSpline(true);
    fScintillatorMPT->AddProperty("RINDEX",PhotonEnergy,rIndex,nEntries);
    fScintillatorMPT->AddProperty("ABSLENGTH",PhotonEnergy,Absorption,nEntries)->SetSpline(true);

    fScintillatorMPT->AddConstProperty("SCINTILLATIONYIELD",60./keV);
    fScintillatorMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
    fScintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 18.0*ns);
    fScintillatorMPT->AddConstProperty("YIELDRATIO",1.0);

// To be add:
//    fScintillatorMPT->AddProperty("WLSABSLENGTH",PhotonEnergy,AbsFiber,nEntries);
//    fScintillatorMPT->AddProperty("WLSCOMPONENT",PhotonEnergy,EmissionFiber,nEntries);


// Reflector optical properties

    G4double rIndex_ref[nEntries] =       {1.35, 1.35, 1.35, 1.35,
                                    1.35, 1.35, 1.35, 1.35, 1.35,
                                    1.35, 1.35, 1.35, 1.35, 1.35,
                                    1.35, 1.35, 1.35, 1.35, 1.35,
                                    1.35, 1.35, 1.35, 1.35, 1.35, 1.35};

    fReflectorMPT->AddProperty("RINDEX",PhotonEnergy,rIndex_ref,nEntries)->SetSpline(true);

// SiPM optical properties

  //  G4double rIndex_SiPM[nEntries] =       {1.35, 1.35, 1.35, 1.35,
  //                                  1.35, 1.35, 1.35, 1.35, 1.35,
  //                                  1.35, 1.35, 1.35, 1.35, 1.35,
  //                                  1.35, 1.35, 1.35, 1.35, 1.35,
  //                                  1.35, 1.35, 1.35, 1.35, 1.35, 1.35};

    G4double Absorption_SiPM[nEntries] =   {0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                                     0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                                     0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                                     0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                                     0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm};

    fSiPMMPT->AddProperty("RINDEX",PhotonEnergy,rIndex,nEntries);
    fSiPMMPT->AddProperty("ABSLENGTH",PhotonEnergy,Absorption_SiPM,nEntries);

// World optical properties
    G4double rIndex_air[nEntries] =       { 1.0003, 1.0003, 1.0003, 1.0003,
                                            1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
                                            1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
                                            1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
                                            1.0003, 1.0003, 1.0003, 1.0003, 1.0003, 1.0003};

    fWorldMPT->AddConstProperty("ABSLENGTH",1 *m);
    fWorldMPT->AddProperty("RINDEX",PhotonEnergy,rIndex_air,nEntries);

//
// ------------ Generate & Add Material Properties Table ------------
//

  Vacuum   ->SetMaterialPropertiesTable(fWorldMPT);
  CeBr_3->SetMaterialPropertiesTable(fScintillatorMPT);
  C2F4  ->SetMaterialPropertiesTable(fReflectorMPT);
  Si    ->SetMaterialPropertiesTable(fSiPMMPT);

  // ------------- Volumes --------------

  // The experimental Hall
  G4Box* world_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);

  G4LogicalVolume* world_LV
    = new G4LogicalVolume(world_box, Vacuum, "World", 0, 0, 0);

  G4VPhysicalVolume* world_PV
    = new G4PVPlacement(0, G4ThreeVector(), world_LV, "World", 0, false, 0);


  // The Air Tank
  G4Box* airTank_box = new G4Box("Tank",fTank_x,fTank_y,fTank_z);

  G4LogicalVolume* airTank_log
    = new G4LogicalVolume(airTank_box,Vacuum,"Tank",0,0,0);

  fTank
    = new G4PVPlacement(0, G4ThreeVector(), airTank_log, "Tank",
                        world_LV, false, 0);


    if (!isQuadratic) {
        // The Scintillator
// ************** cylindrical **************
  G4Tubs* Scintillator_tube
    = new G4Tubs ("Scintillator", 0*m, fScintillator_r, fScintillator_z, 0*deg, 360*deg);

  Scintillator_log
    = new G4LogicalVolume(Scintillator_tube, CeBr_3,"Scintillator",0,0,0);

            // The SiPM

// ************** cylindrical **************
  G4Tubs* SiPM_tube
    = new G4Tubs ("SiPM", 0, fScintillator_r, fSiPM_z, 0*deg, 360*deg);

  SiPM_log
    = new G4LogicalVolume(SiPM_tube, Si,"SiPM",0,0,0);

        if (hasReflector) {

             // The reflector around ----- cylindrical **************************************************
  G4Tubs* Reflector_tube
    = new G4Tubs ("Reflector", fScintillator_r, fScintillator_r+0.3*cm, fScintillator_z+0.3*cm, 0*deg, 360*deg);

  G4LogicalVolume* Reflector_log
    = new G4LogicalVolume(Reflector_tube, C2F4,"Reflector",0,0,0);

  // visibility of the Reflector
  G4VisAttributes*VisAtt_Reflector_log
    = new G4VisAttributes(true, G4Color(0.7, 0.0, 0.3, 0.8));

  VisAtt_Reflector_log->SetForceSolid(true);
  Reflector_log->SetVisAttributes(VisAtt_Reflector_log);

  fReflector
    = new G4PVPlacement(0, G4ThreeVector(0,0,0), Reflector_log, "Reflector",
                        airTank_log, false, 0);


    // The reflector infront   *******************************************************
  G4Tubs* Reflector_tube_2
    = new G4Tubs ("Reflector_2", 0, fScintillator_r, fSiPM_z, 0*deg, 360*deg);

  G4LogicalVolume* Reflector_log_2
    = new G4LogicalVolume(Reflector_tube_2, C2F4,"Reflector_2",0,0,0);

  // visibility of the Reflector_2
  G4VisAttributes*VisAtt_Reflector_log_2
    = new G4VisAttributes(true, G4Color(0.7, 0.0, 0.4, 0.8));

  VisAtt_Reflector_log_2->SetForceSolid(true);
  Reflector_log_2->SetVisAttributes(VisAtt_Reflector_log_2);

  fReflector_2
    = new G4PVPlacement(0, G4ThreeVector(0,0,-fScintillator_z-fSiPM_z), Reflector_log_2, "Reflector_2",
                        airTank_log, false, 0);
        }
    }

    if (isQuadratic) {
// ************** quadratic scintillator **************
  G4Box* Scintillator_box
    = new G4Box ("Scintillator", fScintillator_a, fScintillator_a, fScintillator_z);

  Scintillator_log
    = new G4LogicalVolume(Scintillator_box, CeBr_3,"Scintillator",0,0,0);

  // ************** quadratic SiPM **************
  G4Box* SiPM_box
    = new G4Box ("SiPM", fScintillator_a, fScintillator_a, fSiPM_z);

  SiPM_log
    = new G4LogicalVolume(SiPM_box, Si,"SiPM",0,0,0);

        if (hasReflector) {
       // The reflector around -----quadratic**  **************************************************+
    G4Box* Reflector1_box
    = new G4Box ("Reflector1", 0.15*cm, fScintillator_a+0.3*cm, fScintillator_z+0.3*cm);

    G4LogicalVolume* Reflector1_log
    = new G4LogicalVolume(Reflector1_box, C2F4,"Reflector",0,0,0);


    G4Box* Reflector2_box
    = new G4Box ("Reflector2", fScintillator_a, 0.15*cm, fScintillator_z+0.3*cm);

    G4LogicalVolume* Reflector2_log
    = new G4LogicalVolume(Reflector2_box, C2F4,"Reflector",0,0,0);


    G4Box* Reflector3_box
    = new G4Box ("Reflector3", fScintillator_a, fScintillator_a, 0.15*cm);

    G4LogicalVolume* Reflector3_log
    = new G4LogicalVolume(Reflector3_box, C2F4,"Reflector",0,0,0);


  // visibility of the quadratic Reflector
    G4VisAttributes*VisAtt_Reflector_log_123
    = new G4VisAttributes(true, G4Color(0.7, 0.0, 0.3, 0.8));
    VisAtt_Reflector_log_123->SetForceSolid(true);

    Reflector1_log->SetVisAttributes(VisAtt_Reflector_log_123);
    Reflector2_log->SetVisAttributes(VisAtt_Reflector_log_123);
    Reflector3_log->SetVisAttributes(VisAtt_Reflector_log_123);



			//Physical volume --- quadratic ****
    fReflector11
    = new G4PVPlacement(0, G4ThreeVector(fScintillator_a + 0.15*cm,0,0), Reflector1_log, "Reflector11",
                        airTank_log, false, 0);

    fReflector12
    = new G4PVPlacement(0, G4ThreeVector(-fScintillator_a - 0.15*cm,0,0), Reflector1_log, "Reflector12",
                        airTank_log, false, 0);

    fReflector21
    = new G4PVPlacement(0, G4ThreeVector(0,fScintillator_a + 0.15*cm,0), Reflector2_log, "Reflector21",
                        airTank_log, false, 0);

    fReflector22
    = new G4PVPlacement(0, G4ThreeVector(0,-fScintillator_a - 0.15*cm,0), Reflector2_log, "Reflector22",
                        airTank_log, false, 0);

    fReflector3
    = new G4PVPlacement(0, G4ThreeVector(0,0,-fScintillator_z - 0.15*cm), Reflector3_log, "Reflector3",
                        airTank_log, false, 0);
        }
    }

    //********* visibility **********

    // visibility of the Scintillator
  G4VisAttributes*VisAtt_Scintillator_log
    = new G4VisAttributes(true, G4Color(0.4, 0.6, 0.1, 1));

  VisAtt_Scintillator_log->SetForceSolid(true);
  Scintillator_log->SetVisAttributes(VisAtt_Scintillator_log);

  G4VPhysicalVolume* fScintillator
    = new G4PVPlacement(0, G4ThreeVector(0,0,0), Scintillator_log, "Scintillator",
                        airTank_log, false, 0);

  // visibility of the SiPM
  G4VisAttributes*VisAtt_SiPM_log
    = new G4VisAttributes(true, G4Color(0.2, 0.0, 0.8, 0.8));

  VisAtt_SiPM_log->SetForceSolid(true);
  SiPM_log->SetVisAttributes(VisAtt_SiPM_log);

  G4VPhysicalVolume* fSiPM
    = new G4PVPlacement(0, G4ThreeVector(0,0,fScintillator_z+fSiPM_z), SiPM_log, "SiPM",
                        airTank_log, false, 0);

  // ------------- Surface --------------

    if (hasReflector) {
  // Scintillator - Reflector surface  ***************************************************************

    G4OpticalSurface* Sc_Ref = new G4OpticalSurface("Sc_Ref");

    G4LogicalBorderSurface* Sc_Ref_log; //cylindrical
    G4LogicalBorderSurface* Sc_Ref_log_2; //cylindrical

    G4LogicalBorderSurface* Sc_Ref_log_11; //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_12;  //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_21; //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_22;  //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_3;  //quadratic

    if (isQuadratic) {
        Sc_Ref_log_11      = new G4LogicalBorderSurface("Sc_Ref11", fScintillator, fReflector11, Sc_Ref);
        Sc_Ref_log_12    = new G4LogicalBorderSurface("Sc_Ref12", fScintillator, fReflector12, Sc_Ref);
        Sc_Ref_log_21      = new G4LogicalBorderSurface("Sc_Ref21", fScintillator, fReflector21, Sc_Ref);
        Sc_Ref_log_22    = new G4LogicalBorderSurface("Sc_Ref22", fScintillator, fReflector22, Sc_Ref);
        Sc_Ref_log_3    = new G4LogicalBorderSurface("Sc_Ref3", fScintillator, fReflector3, Sc_Ref);

    }
    else {
        Sc_Ref_log      = new G4LogicalBorderSurface("Sc_Ref", fScintillator, fReflector, Sc_Ref);
        Sc_Ref_log_2    = new G4LogicalBorderSurface("Sc_Ref", fScintillator, fReflector_2, Sc_Ref);
    }

    Sc_Ref->SetType(dielectric_dielectric);
    Sc_Ref->SetModel(unified);
    if (Reflector == 2) {
		Sc_Ref->SetFinish(polishedbackpainted);
    }
    else {
		Sc_Ref->SetFinish(groundbackpainted);
    }
    if (Reflector == 3) {
	    Sc_Ref->SetSigmaAlpha(0.6);
    }
    else {
	    Sc_Ref->SetSigmaAlpha(0.0);
    }


    G4double reflectivity[nEntries]=      { 0.97, 0.97, 0.97, 0.97, 0.97,
                                            0.97, 0.97, 0.97, 0.97, 0.97,
                                            0.97, 0.97, 0.97, 0.97, 0.97,
                                            0.97, 0.97, 0.97, 0.97, 0.97,
                                            0.97, 0.97, 0.97, 0.97, 0.97};


    G4double SpSp[nEntries]=       	{1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0};



    if (Reflector == 0) {
	fSc_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex, nEntries);
    }
    else {
	fSc_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex_air, nEntries);
    }
    fSc_RefMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity, nEntries);
    if (Reflector == 3){
	    fSc_RefMPT->AddProperty("LAMBERTIAN", PhotonEnergy, SpSp, nEntries);
    }
    else {
    	fSc_RefMPT->AddProperty("SPECULARSPIKECONSTANT", PhotonEnergy, SpSp, nEntries);
    }
    Sc_Ref->SetMaterialPropertiesTable(fSc_RefMPT);


// Scintillator - SiPM surface   ***************************************************************************

//    G4OpticalSurface* Sc_SiPM = new G4OpticalSurface("Sc_SiPM");

//    G4LogicalBorderSurface* Sc_SiPM_log = new G4LogicalBorderSurface("Sc_SiPM", fScintillator, fSiPM, Sc_SiPM);

//    Sc_SiPM->SetType(dielectric_dielectric);
//    Sc_SiPM->SetModel(unified);
//    Sc_SiPM->SetFinish(polished);



//    G4double reflectivity_SiPM[nEntries]=      { 1.0, 1.0, 1.0, 1.0,
//                                            1.0, 1.0, 1.0, 1.0,
//                                            1.0, 1.0, 1.0, 1.0,
//                                            1.0, 1.0, 1.0, 1.0,
//                                            1.0, 1.0, 1.0, 1.0, 1.0};


//    fSc_SiPMMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity_SiPM, nEntries);
//   fSc_SiPMMPT->AddProperty("RINDEX", PhotonEnergy, rIndex, nEntries);

//    Sc_SiPM->SetMaterialPropertiesTable(fSc_SiPMMPT);


// SiPM - Reflector surface             *************************************************+

//    G4OpticalSurface* Si_Ref = new G4OpticalSurface("Si_Ref");

//    G4LogicalBorderSurface* Si_Ref_log      = new G4LogicalBorderSurface("Si_Ref", fSiPM, fReflector, Si_Ref);
//    G4LogicalBorderSurface* Si_Ref_log_11      = new G4LogicalBorderSurface("Si_Ref", fSiPM, fReflector11, Si_Ref);
//    G4LogicalBorderSurface* Si_Ref_log_12      = new G4LogicalBorderSurface("Si_Ref", fSiPM, fReflector12, Si_Ref);
//    G4LogicalBorderSurface* Si_Ref_log_21      = new G4LogicalBorderSurface("Si_Ref", fSiPM, fReflector21, Si_Ref);
//    G4LogicalBorderSurface* Si_Ref_log_22      = new G4LogicalBorderSurface("Si_Ref", fSiPM, fReflector22, Si_Ref);

//    Si_Ref->SetType(dielectric_dielectric);
//    Si_Ref->SetModel(unified);
//    Si_Ref->SetFinish(groundbackpainted);
//    Si_Ref->SetSigmaAlpha(0.0);

//    fSi_RefMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity, nEntries);
//    fSi_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex, nEntries);
//    fSi_RefMPT->AddProperty("LAMBERTIAN", PhotonEnergy, SpSp, nEntries);
//    Si_Ref->SetMaterialPropertiesTable(fSi_RefMPT);

    }






    return world_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfaceSigmaAlpha(G4double v) {
  fSurface->SetSigmaAlpha(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface sigma alpha set to: " << fSurface->GetSigmaAlpha()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddBoxMPV(const char* c,
                                     G4MaterialPropertyVector* mpv) {
  mpv->SetSpline(true);
  fBoxMPT->AddProperty(c, mpv);
  G4cout << "The MPT for the box is now: " << G4endl;
  fBoxMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPV(const char* c,
                                       G4MaterialPropertyVector* mpv) {
  mpv->SetSpline(true);
  fWorldMPT->AddProperty(c, mpv);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPV(const char* c,
                                         G4MaterialPropertyVector* mpv) {
  mpv->SetSpline(true);
  fSurfaceMPT->AddProperty(c, mpv);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddBoxMPCV(const char* c, G4double v) {
  fBoxMPT->AddConstProperty(c, v);
  G4cout << "The MPT for the box is now: " << G4endl;
  fBoxMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPCV(const char* c, G4double v) {
  fWorldMPT->AddConstProperty(c, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;

}
