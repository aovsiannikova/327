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

#include <ActionInitialization.hh>
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "DetectorMaterials.hh"

#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"
#include "G4Para.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4int detectorForm, G4bool hasReflector, G4bool hasAlu,
                                           G4double ScintillatorThickness, G4int ReflectorType,
                                           G4double ScintillatorSideLength, G4double ScintillatorRadius,
                                           G4int ScintillatorMaterial, G4double ScintillatorWidth,
                                           G4int NumberOfSiPM
)
 : G4VUserDetectorConstruction(),
   fScintillatorThickness(ScintillatorThickness),
   fScintillatorLength(ScintillatorSideLength),
   fScintillatorWidth(ScintillatorWidth),

   fScintillatorRadius(ScintillatorRadius),
   fScintillatorMaterial(ScintillatorMaterial),
   fDetectorSection(detectorForm),

   fHasReflector(hasReflector),
   fHasAlu(hasAlu),
   fNumberOfSiPM(NumberOfSiPM),
   fDetectorMessenger(nullptr),
   fReflectorType(ReflectorType)

{
  fSiPM_z = 1.5*mm;
  fRefl_z = 1.5*mm;
  fAlu_thickness = 0. * mm;

  fSurface = nullptr;
  fScintillator = nullptr;

  fWorldMPT           = new G4MaterialPropertiesTable();
  fScintillatorMPT    = new G4MaterialPropertiesTable();
  fReflectorMPT       = new G4MaterialPropertiesTable();
  fSiPMMPT            = new G4MaterialPropertiesTable();
  fSurfaceMPT         = new G4MaterialPropertiesTable();
  fSc_RefMPT          = new G4MaterialPropertiesTable();

  fDetectorMessenger  = new DetectorMessenger(this);
  DefineMaterials();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (fDetectorSection == 1) {
    fExpHall_x = fExpHall_y = fExpHall_z =
            (std::max({fScintillatorThickness*2, fScintillatorLength*2, fScintillatorWidth*2}) + 10.) * mm;
  }
  else if (fDetectorSection == 2) {
    fExpHall_x = fExpHall_y = fExpHall_z =
            (std::max({fScintillatorThickness*2, fScintillatorLength*2, fScintillatorLength*2}) + 10.) * mm;
  }
  else {
    fExpHall_x = fExpHall_y = fExpHall_z =
            (std::max({fScintillatorThickness, fScintillatorRadius}) + 10.) * mm;
  }

  // ------------- Add optical properties -------------
  DefineOpticalProperties(fScintillatorMaterial);

  // ------------- Add Volumes --------------

  // The experimental Hall
  auto* world_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);

  auto* world_LV
    = new G4LogicalVolume(world_box, Vacuum, "World", nullptr, nullptr, nullptr);

  G4VPhysicalVolume* world_PV
    = new G4PVPlacement(nullptr, G4ThreeVector(), world_LV, "World",
                        nullptr, false, 0);


    // visibility
  auto* VisAtt_Reflector_log
          = new G4VisAttributes(true, G4Color(0.7, 0.0, 0.3, 1.0));
  auto* VisAtt_Alu_log
          = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.4, 1.0));
  auto* VisAtt_Scintillator_log
          = new G4VisAttributes(true, G4Color(0.4, 0.6, 0.1, 1));
  auto* VisAtt_SiPM_log
          = new G4VisAttributes(true, G4Color(0.2, 0.0, 0.8, 1.0));

  VisAtt_SiPM_log->SetForceSolid(true);
  VisAtt_Reflector_log->SetForceSolid(true);
  VisAtt_Alu_log->SetForceSolid(true);
  VisAtt_Scintillator_log->SetForceSolid(true);

  auto* Sc_Ref = new G4OpticalSurface("Sc_Ref");

// ************** quadratic  **************
  if (fDetectorSection == 1) {

      // ************** quadratic scintillator **************
    AddQScintillator(VisAtt_Scintillator_log, world_LV, fScintillatorLength, fScintillatorWidth, fScintillatorThickness, fScintillatorMaterial);

    if (fHasReflector) {
      if (fHasAlu)           {
        fAlu_thickness = 0.25 * mm;
        fRefl_z=0.5*mm;
        fSiPM_z = 0.75*mm;
                 // ************** quadratic SiPM **************

        AddQSiPM(VisAtt_SiPM_log, world_LV, fNumberOfSiPM, fScintillatorWidth, fScintillatorLength, fScintillatorThickness, fSiPM_z);
        // ************ quadratic reflector
        AddQReflector(VisAtt_Reflector_log, world_LV, fNumberOfSiPM, fScintillatorLength, fScintillatorWidth, fSiPM_z, fRefl_z);

        // ************ quadratic Aluminium housing
        AddQHousing(VisAtt_Alu_log, world_LV, fNumberOfSiPM, fScintillatorLength + 2 * fRefl_z,
                    fScintillatorWidth+2* fRefl_z, fScintillatorThickness + fRefl_z);
      }

      else {
        // ************** quadratic SiPM **************
        AddQSiPM(VisAtt_SiPM_log, world_LV, fNumberOfSiPM, fScintillatorWidth, fScintillatorLength, fScintillatorThickness, fSiPM_z);

        // ************ quadratic reflector
        AddQReflector(VisAtt_Reflector_log, world_LV, fNumberOfSiPM, fScintillatorLength, fScintillatorWidth, fSiPM_z, fRefl_z);
      }
      // ------------- Surface --------------
      Sc_RefLogYZ1      = new G4LogicalBorderSurface("Sc_RefYZ1", fScintillator, fReflectorYZ1, Sc_Ref);
      Sc_RefLogYZ2    = new G4LogicalBorderSurface("Sc_RefYZ2", fScintillator, fReflectorYZ2, Sc_Ref);
      Sc_Ref_log_21      = new G4LogicalBorderSurface("Sc_Ref21", fScintillator, fReflectorXZ1, Sc_Ref);
      Sc_Ref_log_22    = new G4LogicalBorderSurface("Sc_Ref22", fScintillator, fReflectorXZ2, Sc_Ref);
      Sc_Ref_log_3    = new G4LogicalBorderSurface("Sc_Ref3", fScintillator, fReflectorXY1, Sc_Ref);
      ChooseReflector(fReflectorType, Sc_Ref, fScintillatorMaterial);
    }

    else {
      // ************** quadratic SiPM **************
      AddQSiPM(VisAtt_SiPM_log, world_LV, fNumberOfSiPM, fScintillatorWidth, fScintillatorLength, fScintillatorThickness, fSiPM_z);

      if (fHasAlu) {
        fAlu_thickness = 0.25 * mm;
        fRefl_z=0.*mm;
        // ************ quadratic Aluminium housing
        AddQHousing(VisAtt_Alu_log, world_LV, fNumberOfSiPM, fScintillatorLength, fScintillatorWidth, fScintillatorThickness);
      }
    }
    }
  // ************** cylindrical
  else if (fDetectorSection == 0) {
    // ************** cylindrical Scintillator **************
    AddCScintillator(VisAtt_Scintillator_log, world_LV, fScintillatorRadius, fScintillatorThickness, fScintillatorMaterial);

    if (fHasReflector) {
      if (fHasAlu) {
        fRefl_z=0.5*mm;
        fAlu_thickness = 0.25 * mm;
        fSiPM_z = 0.75*mm;

        // ************** cylindrical SiPM **************
        AddCSiPM(VisAtt_SiPM_log, world_LV, fScintillatorRadius, fSiPM_z);

        // The Alu housing infront and around  *******************************************************
        AddCHousing(VisAtt_Alu_log, world_LV, fScintillatorRadius+fRefl_z, fSiPM_z);

        AddCReflector(VisAtt_Reflector_log, world_LV, fScintillatorRadius, fSiPM_z, fRefl_z);


      }

      else {
        // ************** cylindrical SiPM **************
        AddCSiPM(VisAtt_SiPM_log, world_LV, fScintillatorRadius, fSiPM_z);
        AddCReflector(VisAtt_Reflector_log, world_LV, fScintillatorRadius, fSiPM_z, fRefl_z);
        }
      // ------------- Surface --------------
      Sc_Ref_log      = new G4LogicalBorderSurface("Sc_Ref", fScintillator, fReflector, Sc_Ref);
      Sc_Ref_log_2    = new G4LogicalBorderSurface("Sc_Ref", fScintillator, fReflector_2, Sc_Ref);
      ChooseReflector(fReflectorType, Sc_Ref, fScintillatorMaterial);
    }
    else {
      AddCSiPM(VisAtt_SiPM_log, world_LV, fScintillatorRadius, fSiPM_z);

      if (fHasAlu) {
        fAlu_thickness = 0.25 * mm;

        // The Alu infront   *******************************************************
        AddCHousing(VisAtt_Alu_log, world_LV, fScintillatorRadius, fSiPM_z);
      }
    }
  }

  // ************** hexagonal Scintillator **************
  else if (fDetectorSection == 2){
    x0=fScintillatorLength;
    y0=fScintillatorLength;
    fSiPM_z = 0.75*mm;

    if (fHasReflector) {
      if (fHasAlu){
        fAlu_thickness = 0.25 * mm;
        fRefl_z=0.5*mm;
      }
      else{
        fRefl_z=fSiPM_z;
      }
    }
    else {
      fRefl_z = 0.*mm;
    }

    x1=sqrt(3)*fSiPM_z*sqrt(3)/2;
    x2=sqrt(3) * fRefl_z;
    x4=sqrt(3)*x0;
    x3=x4+fSiPM_z;
    x5=x4+fRefl_z;

    y1=2*x0+2*fSiPM_z/sqrt(3);
    y2=y0+fSiPM_z*sqrt(3)/3;
    y3=y0+fRefl_z;
    y4=2*x0+2*fRefl_z/sqrt(3);
    y6=y0+fRefl_z*sqrt(3)/3;

    if (fScintillatorMaterial==0) {
        AddHScintillator(VisAtt_Scintillator_log, world_LV, fScintillatorLength, fScintillatorThickness, CeBr_3);
      }
    else if (fScintillatorMaterial==1) {
        AddHScintillator(VisAtt_Scintillator_log, world_LV, fScintillatorLength, fScintillatorThickness, CsI_Tl);
      }
    else if (fScintillatorMaterial==2) {
        AddHScintillator(VisAtt_Scintillator_log, world_LV, fScintillatorLength, fScintillatorThickness, NaI_Tl);
      }
    AddHSiPM(VisAtt_SiPM_log, world_LV, fScintillatorLength, fScintillatorThickness, fSiPM_z);

    if (fHasReflector) {
      // ************ reflector *
      AddHReflector(VisAtt_Reflector_log, world_LV, fScintillatorLength, fScintillatorThickness, fRefl_z);

      // ------------- Surfaces --------------
      Sc_RefLogTr11      = new G4LogicalBorderSurface("Sc_RefTr11", fScintillator_B, fReflectorTrap11, Sc_Ref);
      Sc_RefLogTr12      = new G4LogicalBorderSurface("Sc_RefTr12", fScintillator_B, fReflectorTrap12, Sc_Ref);
      Sc_RefLogTr21      = new G4LogicalBorderSurface("Sc_RefTr21", fScintillator_C, fReflectorTrap21, Sc_Ref);
      Sc_RefLogTr22      = new G4LogicalBorderSurface("Sc_RefTr22", fScintillator_C, fReflectorTrap22, Sc_Ref);

      Sc_RefLogBox1    = new G4LogicalBorderSurface("Sc_RefBox1", fScintillator_A, fReflectorBox1, Sc_Ref);
      Sc_RefLogBox2    = new G4LogicalBorderSurface("Sc_RefBox2", fScintillator_A, fReflectorBox2, Sc_Ref);

      Sc_RefLogSide1      = new G4LogicalBorderSurface("Sc_RefSide1", fScintillator_B, fReflectorS1, Sc_Ref);
      Sc_RefLogSide2      = new G4LogicalBorderSurface("Sc_RefSide2", fScintillator_C, fReflectorS2, Sc_Ref);

      ChooseReflector(fReflectorType, Sc_Ref, fScintillatorMaterial);
    }
    if (fHasAlu){
      // ************ hexagonal Aluminium housing
      AddHHousing(VisAtt_Alu_log, world_LV, fScintillatorLength, fAlu_thickness, fScintillatorThickness);
    }

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
void DetectorConstruction::SetCrystalThickness(G4double v) {
  fScintillatorThickness=v/2;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Crystal thickness set to: " << v << " mm"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetQuadraticCrystalSideLength(G4double v) {
  fScintillatorLength= v / 2;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Crystal length set to: " << v << " mm"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetQuadraticCrystalWidth(G4double v) {
  fScintillatorWidth= v / 2;
  G4RunManager::GetRunManager()->ReinitializeGeometry();

  G4cout << "Crystal width set to: " << v << " mm"
         << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCylindricalCrystalRadius(G4double v) {
  fScintillatorRadius=v/2;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Crystal radius set to: " << v << " mm"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetReflectorType(G4int v) {
  fReflectorType=v;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Reflector type set to: " << fReflectorType
         << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetectorSection(G4int v) {
  G4String section;
  if (v>2) {
    if (fDetectorSection==0) {section ="Cylindrical";}
    else if (fDetectorSection==1) {section ="Quadratic";}
    else if (fDetectorSection==2) {section ="Hexagonal";}
    G4cout << "Unknown section can not be set. " << section << " section still will be used." << G4endl;
  }
  else {
    fDetectorSection = v;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    if (v==0) {section ="cylindrical";}
    else if (v==1) {section ="quadratic";}
    else if (v==2) {section ="hexagonal";}
    G4cout << "Detector is now "<< section <<"."<< G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetectorHasReflector(G4bool v) {
  fHasReflector = v;

  G4RunManager::GetRunManager()->ReinitializeGeometry();
  if (v==1) {
    G4cout << "Reflector has been added"  << G4endl;
  }
  else {
    G4cout << "Detector has no reflector"  << G4endl;
  }
}

void DetectorConstruction::SetScintillatorMaterial(G4int v){
  G4String material;

  if (v>2) {
    if (fScintillatorMaterial==0) {material ="CeBr3";}
    else if (fScintillatorMaterial==1) {material ="CsI(Tl)";}
    else if (fScintillatorMaterial==2) {material ="NaI(Tl)";}

    G4cout << "Unknown material. " << fScintillatorMaterial << " still will be used." << G4endl;
  }
  else {
    fScintillatorMaterial = v;
    if (v==0) {material ="CeBr3";}
    else if (v==1) {material ="CsI(Tl)";}
    else if (v==2) {material ="NaI(Tl)";}
    G4cout << material << " is used now in scintillator." << G4endl;
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetectorHasAlu(G4bool v) {
  fHasAlu = v;

  G4RunManager::GetRunManager()->ReinitializeGeometry();
  if (v==1) {
    G4cout << "Alu housing  has been added"  << G4endl;
  }
  else {
    G4cout << "Detector has no housing"  << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ChangeNumberOfSiPM(G4int v) {
  fNumberOfSiPM = v;

  G4RunManager::GetRunManager()->ReinitializeGeometry();
  if (v==2) {
    G4cout << "Detector has 2 SiPM (XY, YZ2)"  << G4endl;
  }
  else if (v==1) {
    G4cout << "Detector has one SiPM(XY)"  << G4endl;
  }
  else if (v==3) {
    G4cout << "Detector has 3 SiPM (XY, YZ2, YZ1)"  << G4endl;
  }
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

void DetectorConstruction::DefineMaterials() {
  // ------------- Create Materials -------------
  // CeBr3 (the scintillator material)
  CeBr_3 = new G4Material("CeBr_3", 5.2*g/cm3, 2);
  CeBr_3->AddElement(elCe, 0.36888);
  CeBr_3->AddElement(elBr, 0.63112);

  // CsI (the scintillator material 1)
  CsI = new G4Material("CsI", 4.51*g/cm3, 2);
  CsI->AddElement(elCs, 0.511547);
  CsI->AddElement(elI, 0.488453);

  //CsI(Tl)
  CsI_Tl = new G4Material("CsI_Tl", 4.51*g/cm3, 2);
  CsI_Tl->AddMaterial(CsI,99.6*perCent);
  CsI_Tl->AddElement(elTl,0.4*perCent);

  // NaI (the scintillator material 2)
  NaI = new G4Material("NaI", 3.67*g/cm3, 2);
  NaI->AddElement(elNa, 0.153436);
  NaI->AddElement(elI, 0.846564);

  //NaI(Tl)
  NaI_Tl = new G4Material("NaI_Tl", 3.67*g/cm3, 2);
  NaI_Tl->AddMaterial(NaI,99.6*perCent);
  NaI_Tl->AddElement(elTl,0.4*perCent);

  // vacuum - tank, world
  Vacuum = new G4Material("Vacuum", z = 1., a = 1.008*g/mole, 1.e-25*g/cm3);

  // Teflon (the reflector material)
  C2F4 = new G4Material("C2F4", 2.2*g/cm3, 2);
  C2F4->AddElement(elC, 0.24);
  C2F4->AddElement(elF, 0.76);
}

void DetectorConstruction::DefineOpticalProperties(G4int ScintillatorMaterial) {
  // Add scintillator optical properties
  if (ScintillatorMaterial==0) {
    // Add CeBr3 optical properties
    fScintillatorMPT->AddProperty("FASTCOMPONENT", PhotonEnergyCeBr3, FastCompCeBr3, nEntries)->SetSpline(true);
    fScintillatorMPT->AddProperty("RINDEX", PhotonEnergyCeBr3, rIndexCeBr3, nEntries);
    fScintillatorMPT->AddProperty("ABSLENGTH", PhotonEnergyCeBr3, AbsorptionCeBr3, nEntries)->SetSpline(true);

    fScintillatorMPT->AddConstProperty("SCINTILLATIONYIELD", 60. / keV);
    fScintillatorMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fScintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 18.0 * ns);
    fScintillatorMPT->AddConstProperty("YIELDRATIO", 1.0);

    // Add reflector optical properties
    fReflectorMPT->AddProperty("RINDEX", PhotonEnergyCeBr3, rIndex_ref, nEntries)->SetSpline(true);

// To be add:
//    fScintillatorMPT->AddProperty("WLSABSLENGTH",PhotonEnergyCeBr3,AbsFiber,nEntries);
//    fScintillatorMPT->AddProperty("WLSCOMPONENT",PhotonEnergyCeBr3,EmissionFiber,nEntries);
  }
  else if (ScintillatorMaterial==1) {
    // Add CsI(Tl) optical properties
    fScintillatorMPT->AddProperty("FASTCOMPONENT", PhotonEnergyCsI, FastCompCsI, nEntries)->SetSpline(true);
    fScintillatorMPT->AddProperty("RINDEX", PhotonEnergyCsI, rIndexCsI, nEntries);
    fScintillatorMPT->AddProperty("ABSLENGTH", PhotonEnergyCsI, AbsorptionCsI, nEntries)->SetSpline(true);

    fScintillatorMPT->AddConstProperty("SCINTILLATIONYIELD", 54 / keV);
    fScintillatorMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fScintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 0.6 * us);
    // fScintillatorMPT->AddConstProperty("SLOWTIMECONSTANT", 3.5 * mus);
    fScintillatorMPT->AddConstProperty("YIELDRATIO", 1.0);

    // Add reflector optical properties
    fReflectorMPT->AddProperty("RINDEX", PhotonEnergyCsI, rIndex_ref, nEntries)->SetSpline(true);

  }
  else if (ScintillatorMaterial==2) {
    // Add NaI(Tl) optical properties
    fScintillatorMPT->AddProperty("FASTCOMPONENT", PhotonEnergyNaI, FastCompNaI, nEntries)->SetSpline(true);
    fScintillatorMPT->AddProperty("RINDEX", PhotonEnergyNaI, rIndexNaI, nEntries);
    fScintillatorMPT->AddProperty("ABSLENGTH", PhotonEnergyNaI, AbsorptionNaI, nEntries)->SetSpline(true);

    fScintillatorMPT->AddConstProperty("SCINTILLATIONYIELD", 38 / keV);
    fScintillatorMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fScintillatorMPT->AddConstProperty("FASTTIMECONSTANT", 250.0 * ns);
    fScintillatorMPT->AddConstProperty("YIELDRATIO", 1.0);

    // Add reflector optical properties
    fReflectorMPT->AddProperty("RINDEX", PhotonEnergyNaI, rIndex_ref, nEntries)->SetSpline(true);

  }


// Add SiPM optical properties
  if (ScintillatorMaterial==0) {
    fSiPMMPT->AddProperty("RINDEX", PhotonEnergyCeBr3, rIndexCeBr3, nEntries);
    fSiPMMPT->AddProperty("ABSLENGTH", PhotonEnergyCeBr3, Absorption_SiPM, nEntries);
  }
  else if (ScintillatorMaterial==1) {
    fSiPMMPT->AddProperty("RINDEX", PhotonEnergyCsI, rIndexCsI, nEntries);
    fSiPMMPT->AddProperty("ABSLENGTH", PhotonEnergyCsI, Absorption_SiPM, nEntries);
  }
  else if (ScintillatorMaterial==2) {
    fSiPMMPT->AddProperty("RINDEX", PhotonEnergyNaI, rIndexNaI, nEntries);
    fSiPMMPT->AddProperty("ABSLENGTH", PhotonEnergyNaI, Absorption_SiPM, nEntries);
  }


// Add World optical properties
  fWorldMPT->AddConstProperty("ABSLENGTH",1 *m);
  if (ScintillatorMaterial==0) {
    fWorldMPT->AddProperty("RINDEX", PhotonEnergyCeBr3, rIndex_air, nEntries);
  }
  else if (ScintillatorMaterial==1) {
    fWorldMPT->AddProperty("RINDEX", PhotonEnergyCsI, rIndex_air, nEntries);
  }
  else if (ScintillatorMaterial==2) {
    fWorldMPT->AddProperty("RINDEX", PhotonEnergyNaI, rIndex_air, nEntries);
  }

  // ------------ Generate & Add Material Properties Table ------------
  Vacuum   ->SetMaterialPropertiesTable(fWorldMPT);
  if (ScintillatorMaterial==0) {
    CeBr_3->SetMaterialPropertiesTable(fScintillatorMPT);
  }
  else if (ScintillatorMaterial==1) {
    CsI_Tl->SetMaterialPropertiesTable(fScintillatorMPT);
  }
  else if (ScintillatorMaterial==2) {
    NaI_Tl->SetMaterialPropertiesTable(fScintillatorMPT);
  }

  C2F4  ->SetMaterialPropertiesTable(fReflectorMPT);
  Si    ->SetMaterialPropertiesTable(fSiPMMPT);
}


void DetectorConstruction::AddQHousing(G4VisAttributes* VisAtt_Alu_log, G4LogicalVolume* world_LV, G4int NumberOfSiPM,
                                       double AluLength, double AluWidth, double Thickness) {
  auto *AluBoxXY
          = new G4Box("AluBoxXY", AluWidth, AluLength, fAlu_thickness);

  auto *AluBoxXZ
          = new G4Box("AluBoxXY", AluWidth, fAlu_thickness, Thickness+fAlu_thickness+fSiPM_z);

  auto *AluBoxYZ
          = new G4Box("AluBoxXY", fAlu_thickness, AluLength + 2 * fAlu_thickness,
                      Thickness+fAlu_thickness+fSiPM_z);

  AluBoxLogXY
          = new G4LogicalVolume(AluBoxXY, Al, "Alu", nullptr, nullptr, nullptr);
  AluBoxLogXZ
          = new G4LogicalVolume(AluBoxXZ, Al, "Alu", nullptr, nullptr, nullptr);
  AluBoxLogYZ
          = new G4LogicalVolume(AluBoxYZ, Al, "Alu", nullptr, nullptr, nullptr);

// visibility of the quadratic Alu housing
  AluBoxLogXY->SetVisAttributes(VisAtt_Alu_log);
  AluBoxLogXZ->SetVisAttributes(VisAtt_Alu_log);
  AluBoxLogYZ->SetVisAttributes(VisAtt_Alu_log);

  fAluBoxXY1
          = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -fAlu_thickness - Thickness-fRefl_z), AluBoxLogXY,
                              "AluBoxXY1", world_LV, false, 0);
  fAluBoxXZ1
          = new G4PVPlacement(nullptr, G4ThreeVector(0, AluLength+fAlu_thickness, -Thickness+fScintillatorThickness+fSiPM_z-fAlu_thickness), AluBoxLogXZ,
                              "AluBoxXZ1", world_LV, false, 0);
  fAluBoxXZ2
          = new G4PVPlacement(nullptr, G4ThreeVector(0, -AluLength-fAlu_thickness, -Thickness+fScintillatorThickness+fSiPM_z-fAlu_thickness), AluBoxLogXZ,
                              "AluBoxXZ2", world_LV, false, 0);
  if (NumberOfSiPM<3) {
    fAluBoxYZ1
            = new G4PVPlacement(nullptr, G4ThreeVector(AluWidth, 0, -Thickness + fScintillatorThickness + fSiPM_z - fAlu_thickness),
                                AluBoxLogYZ, "AluBoxYZ1", world_LV, false, 0);
  }
  if (NumberOfSiPM<2) {
    fAluBoxYZ2
            = new G4PVPlacement(0, G4ThreeVector(-AluWidth, 0, -Thickness + fScintillatorThickness + fSiPM_z - fAlu_thickness),
                                AluBoxLogYZ, "AluBoxYZ2", world_LV, false, 0);
  }
}


void DetectorConstruction::AddHHousing(G4VisAttributes* VisAtt_Alu_log, G4LogicalVolume* world_LV,
                                       double ScintillatorLength, double AluThickness, double ScintillatorThickness) {
  G4double z_Alu = fScintillatorThickness+2*fRefl_z+AluThickness;

  auto *AluBox
  = new G4Box("AluBox", x3, y2, AluThickness);

  auto * AluTrap1 = new G4GenericTrap ("AluTrap1", AluThickness, std::vector<G4TwoVector>{{-x3,-y2},
                                {x3,-y2},{0,-y1},{0,-y1}, {-x3,-y2},{x3,-y2},{0,-y1}, {0,-y1}});

  auto * AluTrap2 = new G4GenericTrap ("AluTrap2", AluThickness, std::vector<G4TwoVector>{{-x3,y2},
                                   {0,y1},{0,y1},{x3,y2}, {-x3,y2},{0,y1},{0,y1}, {x3,y2}});

  auto * AluSideTrap1 = new G4GenericTrap ("AluSideTrap1", z_Alu+AluThickness,
                                           std::vector<G4TwoVector>{{x3,-y2},
                                 {0,-y1},{0,-y4},{x5,-y6}, {x3,-y2},{0,-y1},{0,-y4}, {x5,-y6}});

  auto * AluSideTrap2 = new G4GenericTrap ("AluSideTrap2", z_Alu+AluThickness,
                                           std::vector<G4TwoVector>{{x3,y2},
                                             {x5,y6},{0,y4},{0,y1}, {x3,y2},{x5,y6},{0,y4}, {0,y1}});


  AluBoxLog
  = new G4LogicalVolume(AluBox, Al, "AluBoxLog", nullptr, nullptr, nullptr);
  AluTrapLog1
  = new G4LogicalVolume(AluTrap1, Al, "AluTrapLog1", nullptr, nullptr, nullptr);
  AluTrapLog2
  = new G4LogicalVolume(AluTrap2, Al, "AluTrapLog2", nullptr, nullptr, nullptr);
  AluSideTrapLog1
  = new G4LogicalVolume(AluSideTrap1, Al, "AluSideTrapLog1", nullptr, nullptr, nullptr);
  AluSideTrapLog2
  = new G4LogicalVolume(AluSideTrap2, Al, "AluSideTrapLog2", nullptr, nullptr, nullptr);

  // visibility of the quadratic Alu housing
  AluBoxLog->SetVisAttributes(VisAtt_Alu_log);
  AluTrapLog1->SetVisAttributes(VisAtt_Alu_log);
  AluTrapLog2->SetVisAttributes(VisAtt_Alu_log);
  AluSideTrapLog1->SetVisAttributes(VisAtt_Alu_log);
  AluSideTrapLog2->SetVisAttributes(VisAtt_Alu_log);


  fAluBox1
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_Alu), AluBoxLog,
                      "AluBox1", world_LV, false, 0);
  fAluBox2
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -z_Alu), AluBoxLog,
                      "AluBox2", world_LV, false, 0);


  fAluTrap11
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_Alu), AluTrapLog1, "AluTrap11", world_LV, false, 0);
  fAluTrap12
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -z_Alu), AluTrapLog1, "AluTrap12", world_LV, false, 0);
  fAluTrap21
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_Alu), AluTrapLog2, "AluTrap21", world_LV, false, 0);
  fAluTrap22
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -z_Alu), AluTrapLog2, "AluTrap22", world_LV, false, 0);

  fAluSideTrap1
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), AluSideTrapLog1, "AluSideTrap1", world_LV, false, 0);
  fAluSideTrap2
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), AluSideTrapLog1, "AluSideTrap2", world_LV, false, 0);

}


void DetectorConstruction::AddHSiPM(G4VisAttributes* VisAtt_SiPM_log, G4LogicalVolume* world_LV,
                                    G4double ScintillatorLength, G4double ScintillatorThickness, G4double SiPMThickness) {

  auto * SiPMSide01 = new G4GenericTrap ("SiPMSide01", ScintillatorThickness, std::vector<G4TwoVector>{{-x3, -y2},
                               {-x3, y2}, {-x4, y0}, {-x4,-y0}, {-x3, -y2}, {-x3, y2}, {-x4, y0}, {-x4,-y0} });
  SiPMLogTrap01   = new G4LogicalVolume(SiPMSide01, Si, "SiPM_Log_tr01", nullptr, nullptr, 0);
  SiPMLogTrap01->SetVisAttributes(VisAtt_SiPM_log);
  fSiPMTrap01   = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                      SiPMLogTrap01, "SiPMTrap01", world_LV, false, 0);

  auto * SiPMSide02 = new G4GenericTrap ("SiPMSide02", ScintillatorThickness, std::vector<G4TwoVector>{{x3, y2},
                        {x3, -y2}, {x4, -y0}, {x4,y0}, {x3, y2}, {x3, -y2}, {x4, -y0}, {x4,y0} });
  SiPMLogTrap02   = new G4LogicalVolume(SiPMSide02, Si, "SiPM_Log_tr02", nullptr, nullptr, 0);
  SiPMLogTrap02->SetVisAttributes(VisAtt_SiPM_log);
  fSiPMTrap02   = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                    SiPMLogTrap02, "SiPMTrap02", world_LV, false, 0);

  auto * SiPMSide11 = new G4GenericTrap ("SiPMSide11", ScintillatorThickness, std::vector<G4TwoVector>{{0, -y1},
                    {-x3, -y2}, {-x4, -y0}, {0,-2*y0}, {0, -y1}, {-x3, -y2}, {-x4, -y0}, {0,-2*y0} });

  auto * SiPMSide12 = new G4GenericTrap ("SiPMSide12", ScintillatorThickness, std::vector<G4TwoVector>{{0, y1},
                    {0,2*y0},{-x4, y0}, {-x3, y2}, {0, y1},  {0, 2*y0}, {-x4, y0}, {-x3, y2}});
  SiPMLogTrap11   = new G4LogicalVolume(SiPMSide11, Si, "SiPM_Log_tr11", nullptr, nullptr, 0);
  SiPMLogTrap12   = new G4LogicalVolume(SiPMSide12, Si, "SiPM_Log_tr12", nullptr, nullptr, 0);

  SiPMLogTrap11->SetVisAttributes(VisAtt_SiPM_log);
  SiPMLogTrap12->SetVisAttributes(VisAtt_SiPM_log);

  fSiPMTrap11
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                      SiPMLogTrap11, "SiPMTrap11", world_LV, false, 0);
  fSiPMTrap12
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                      SiPMLogTrap12, "SiPMTrap12", world_LV, false, 0);
  }

    void DetectorConstruction::AddQSiPM(G4VisAttributes* VisAtt_SiPM_log, G4LogicalVolume* world_LV, G4int NumberOfSiPM,
                                        G4double ScintXSize, G4double ScintYSize, G4double ScintZSize, G4double SiPMThickness) {

  auto *SiPMBoxXY
  = new G4Box("SiPM", ScintXSize, ScintYSize, SiPMThickness);

  SiPMLogXY
  = new G4LogicalVolume(SiPMBoxXY, Si, "SiPM", nullptr, nullptr, 0);
  SiPMLogXY->SetVisAttributes(VisAtt_SiPM_log);

  // **************** SiPM Placement
  fSiPMXY
  = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fScintillatorThickness + SiPMThickness),
                      SiPMLogXY, "SiPMXY", world_LV, false, 0);

  if (NumberOfSiPM>1){
    auto *SiPMBoxYZ
    = new G4Box("SiPM", SiPMThickness, ScintYSize, ScintZSize+SiPMThickness);

    SiPMLogYZ
    = new G4LogicalVolume(SiPMBoxYZ, Si, "SiPM", nullptr, nullptr, 0);
    // **************** SiPM Placement
    SiPMLogYZ->SetVisAttributes(VisAtt_SiPM_log);

    fSiPMYZ2
    = new G4PVPlacement(nullptr, G4ThreeVector(-ScintXSize-SiPMThickness, 0, SiPMThickness),
                        SiPMLogYZ, "SiPMYZ2", world_LV, false, 0);

  }
  if (NumberOfSiPM>2){
    // **************** SiPM Placement
    fSiPMYZ1
    = new G4PVPlacement(nullptr, G4ThreeVector(ScintXSize+SiPMThickness, 0, SiPMThickness),
                        SiPMLogYZ, "SiPMYZ1", world_LV, false, 0);

  }

}


void DetectorConstruction::AddCSiPM(G4VisAttributes* VisAtt_SiPM_log, G4LogicalVolume* world_LV,
                                      G4double SiPMSize, G4double SiPM_z) {
    G4Tubs *SiPM_tube
          = new G4Tubs("SiPM", 0, SiPMSize, SiPM_z, 0 * deg, 360 * deg);

  SiPMLogXY
          = new G4LogicalVolume(SiPM_tube, Si, "SiPM", 0, 0, 0);

  // **************** SiPM Placement
  SiPMLogXY->SetVisAttributes(VisAtt_SiPM_log);

  fSiPMXY
          = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fScintillatorThickness + SiPM_z),
                              SiPMLogXY, "SiPM", world_LV, false, 0);
}

void DetectorConstruction::AddQReflector(G4VisAttributes* VisAtt_Reflector_log, G4LogicalVolume* world_LV, G4int NumberOfSiPM,
                                         G4double ScintillatorSideLength, G4double ScintillatorSideWidth, G4double SiPM_z, G4double Refl_z) {
 // The reflector around -----quadratic**  **************************************************+
auto * ReflectorBoxYZ
        = new G4Box ("Reflector1", Refl_z, ScintillatorSideLength+2*Refl_z, fScintillatorThickness +2*SiPM_z - fAlu_thickness);
auto * ReflectorLogYZ
        = new G4LogicalVolume(ReflectorBoxYZ, C2F4, "ReflectorYZ", 0, 0, 0);
ReflectorLogYZ->SetVisAttributes(VisAtt_Reflector_log);


  if (NumberOfSiPM<3) {

    fReflectorYZ1
            = new G4PVPlacement(0, G4ThreeVector(ScintillatorSideWidth + Refl_z, 0, +fAlu_thickness), ReflectorLogYZ,
                                "ReflectorYZ1", world_LV, false, 0);
  }
if (NumberOfSiPM<2) {
  fReflectorYZ2
          = new G4PVPlacement(0, G4ThreeVector(-ScintillatorSideWidth - Refl_z, 0, +fAlu_thickness),
                              ReflectorLogYZ, "ReflectorYZ2", world_LV, false, 0);
}
auto * ReflectorBoxXZ
        = new G4Box ("Reflector2", ScintillatorSideWidth,
                     Refl_z, fScintillatorThickness +2*SiPM_z - fAlu_thickness);

auto * ReflectorLogXZ
        = new G4LogicalVolume(ReflectorBoxXZ, C2F4, "ReflectorXZ", nullptr, 0, 0);
ReflectorLogXZ->SetVisAttributes(VisAtt_Reflector_log);

fReflectorXZ1
= new G4PVPlacement(0, G4ThreeVector(0,ScintillatorSideLength + Refl_z,+fAlu_thickness),
                    ReflectorLogXZ, "ReflectorXZ1", world_LV, false, 0);
fReflectorXZ2
= new G4PVPlacement(0, G4ThreeVector(0,-ScintillatorSideLength - Refl_z,+fAlu_thickness),
                    ReflectorLogXZ, "ReflectorXZ2", world_LV, false, 0);

auto * ReflectorBoxXY
        = new G4Box ("Reflector3", ScintillatorSideWidth, ScintillatorSideLength, Refl_z);

auto * ReflectorLogXY
        = new G4LogicalVolume(ReflectorBoxXY, C2F4,"Reflector",0,0,0);
ReflectorLogXY->SetVisAttributes(VisAtt_Reflector_log);

  fReflectorXY1 //front
= new G4PVPlacement(0, G4ThreeVector(0,0,-fScintillatorThickness - Refl_z),
                    ReflectorLogXY, "Reflector3", world_LV, false, 0);
}


void DetectorConstruction::AddHReflector(G4VisAttributes* VisAtt_Reflector_log, G4LogicalVolume* world_LV,
        G4double ScintillatorSideLength, G4double ScintillatorThickness, G4double RefThickness) {
  // The reflectors on the sides **************************************************+
  auto * ReflectorSide1 = new G4GenericTrap ("ReflectorSide1", ScintillatorThickness, std::vector<G4TwoVector>{{0,-y4},
                                    {0,-2*y0},{x4,-y0},{x3,-y2}, {0,-y4},{0,-2*y0},{x4,-y0},{x3,-y2}});
  auto * ReflectorSide2 = new G4GenericTrap ("ReflectorSide2", ScintillatorThickness, std::vector<G4TwoVector>{{0,y4},
                                  {x3,y2},{x4,y0},{0,2*y0}, {0,y4},{x3,y2},{x4,y0},{0,2*y0}});

  auto * ReflectorLogS1 = new G4LogicalVolume(ReflectorSide1, C2F4, "ReflectorLogS1", 0, 0, 0);
  ReflectorLogS1->SetVisAttributes(VisAtt_Reflector_log);

  auto * ReflectorLogS2 = new G4LogicalVolume(ReflectorSide2, C2F4, "ReflectorLogS2", 0, 0, 0);
  ReflectorLogS2->SetVisAttributes(VisAtt_Reflector_log);

  fReflectorS1 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), ReflectorLogS1, "ReflectorS1", world_LV, false, 0);
  fReflectorS2 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), ReflectorLogS2, "ReflectorS2", world_LV, false, 0);


  auto * ReflectorBox = new G4Box ("ReflectorBox", x3, y2,
                                   RefThickness);

  auto * ReflectorLogBox = new G4LogicalVolume(ReflectorBox, C2F4, "ReflectorBox", nullptr, 0, 0);
  ReflectorLogBox->SetVisAttributes(VisAtt_Reflector_log);

  fReflectorBox1  = new G4PVPlacement(0, G4ThreeVector(0,0, ScintillatorThickness+RefThickness),
                      ReflectorLogBox, "ReflectorBox1", world_LV, false, 0);
  fReflectorBox2
  = new G4PVPlacement(0, G4ThreeVector(0,0,-ScintillatorThickness-RefThickness),
                      ReflectorLogBox, "ReflectorBox2", world_LV, false, 0);

  auto * ReflectorTrap1 = new G4GenericTrap ("ReflectorTrap1", RefThickness,std::vector<G4TwoVector>{{0,-y1},
                            {0,-y1},{-x3,-y2},{x3,-y2}, {0,-y1},{0,-y1},{-x3,-y2}, {x3,-y2}});
  auto * ReflectorTrap2 = new G4GenericTrap ("ReflectorTrap2", RefThickness,std::vector<G4TwoVector>{{0,y1},
                             {0,y1},{x3,y2},{-x3,y2}, {0,y1},{0,y1},{x3,y2}, {-x3,y2}});

  auto * ReflectorLogTrap1  = new G4LogicalVolume(ReflectorTrap1, C2F4,"ReflectorLogTrap1",0,0,0);
  ReflectorLogTrap1->SetVisAttributes(VisAtt_Reflector_log);

  auto * ReflectorLogTrap2  = new G4LogicalVolume(ReflectorTrap2, C2F4,"ReflectorLogTrap2",0,0,0);
  ReflectorLogTrap2->SetVisAttributes(VisAtt_Reflector_log);

  fReflectorTrap11 = new G4PVPlacement(0, G4ThreeVector(0,0,RefThickness+ScintillatorThickness),
                      ReflectorLogTrap1, "ReflectorTrap11", world_LV, false, 0);
  fReflectorTrap12 = new G4PVPlacement(0, G4ThreeVector(0,0,-RefThickness-ScintillatorThickness),
                                      ReflectorLogTrap1, "ReflectorTrap12", world_LV, false, 0);

  fReflectorTrap21 = new G4PVPlacement(0, G4ThreeVector(0,0,RefThickness+ScintillatorThickness),
                                      ReflectorLogTrap2, "ReflectorTrap21", world_LV, false, 0);
  fReflectorTrap22 = new G4PVPlacement(0, G4ThreeVector(0,0,-RefThickness-ScintillatorThickness),
                                       ReflectorLogTrap2, "ReflectorTrap22", world_LV, false, 0);

}


void DetectorConstruction::AddCHousing(G4VisAttributes* VisAtt_Alu_log, G4LogicalVolume* world_LV, G4double AluLength, G4double SiPM_z){

  G4Tubs* Alu_tube_2         = new G4Tubs ("Alu_2", 0, AluLength, fAlu_thickness, 0 * deg, 360 * deg);

  G4LogicalVolume* Alu_log_2         = new G4LogicalVolume(Alu_tube_2, Al,"Alu_2",0,0,0);
  Alu_log_2->SetVisAttributes(VisAtt_Alu_log);

  fAlu    = new G4PVPlacement(0, G4ThreeVector(0,0, -fScintillatorThickness -2*fRefl_z - fAlu_thickness),
                            Alu_log_2, "Alu_2", world_LV, false, 0);

  // The Alu housing around   *******************************************************
  G4Tubs* Alu_tube_3
        = new G4Tubs ("Alu_3", AluLength, AluLength + fAlu_thickness, fScintillatorThickness + 2 * SiPM_z, 0 * deg, 360 * deg);

  G4LogicalVolume* Alu_log_3
        = new G4LogicalVolume(Alu_tube_3, Al,"Alu_2",0,0,0);
  Alu_log_3->SetVisAttributes(VisAtt_Alu_log);

  fAlu_3
          = new G4PVPlacement(0, G4ThreeVector(0,0,0), Alu_log_3, "Alu_3",
                              world_LV, false, 0);
}

void DetectorConstruction::AddCReflector(G4VisAttributes* VisAtt_Reflector_log, G4LogicalVolume* world_LV,
                                         G4double ScintillatorRadius, G4double SiPM_z, G4double Refl_z){
  // The reflector around placement **************************************************
  auto* Reflector_tube
          = new G4Tubs ("Reflector", ScintillatorRadius, ScintillatorRadius + Refl_z,
                        fScintillatorThickness + SiPM_z + Refl_z, 0 * deg, 360 * deg);

  G4LogicalVolume* Reflector_log
          = new G4LogicalVolume(Reflector_tube, C2F4,"Reflector",0,0,0);
  Reflector_log->SetVisAttributes(VisAtt_Reflector_log);

  fReflector
          = new G4PVPlacement(0, G4ThreeVector(0, 0, SiPM_z - Refl_z), Reflector_log, "Reflector",
                              world_LV, false, 0);


  // The reflector infront   *******************************************************
  G4Tubs* Reflector_tube_2
          = new G4Tubs ("Reflector_2", 0, ScintillatorRadius, Refl_z, 0 * deg, 360 * deg);

  G4LogicalVolume* Reflector_log_2
          = new G4LogicalVolume(Reflector_tube_2, C2F4,"Reflector_2",0,0,0);
  Reflector_log_2->SetVisAttributes(VisAtt_Reflector_log);

  fReflector_2
          = new G4PVPlacement(0, G4ThreeVector(0,0, -fScintillatorThickness - Refl_z),
                              Reflector_log_2, "Reflector_2", world_LV, false, 0);

}


void DetectorConstruction::AddQScintillator(G4VisAttributes* VisAtt_Scintillator_log, G4LogicalVolume* world_LV,
                                           G4double ScintillatorLength, G4double ScintillatorWidth, G4double ScintillatorThickness,
                                           G4int ScintillatorMaterial) {
  auto *Scintillator_box
          = new G4Box("Scintillator", ScintillatorWidth, ScintillatorLength, ScintillatorThickness);

  if (ScintillatorMaterial == 0) {
    Scintillator_log
            = new G4LogicalVolume(Scintillator_box, CeBr_3, "Scintillator", 0, 0, 0);
  }
  else if (ScintillatorMaterial ==1) {
    Scintillator_log
            = new G4LogicalVolume(Scintillator_box, CsI_Tl, "Scintillator", 0, 0, 0);
  }
  else if (ScintillatorMaterial ==2) {
    Scintillator_log
            = new G4LogicalVolume(Scintillator_box, NaI_Tl, "Scintillator", 0, 0, 0);
  }

  // ************* Placement
  Scintillator_log->SetVisAttributes(VisAtt_Scintillator_log);

  fScintillator
          = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Scintillator_log, "Scintillator",
                              world_LV, false, 0);

}


void DetectorConstruction::AddHScintillator(G4VisAttributes* VisAtt_Scintillator_log, G4LogicalVolume* world_LV,
                                            G4double ScintillatorLength, G4double ScintillatorThickness,
                                            G4Material* ScintillatorMaterial) {
  auto *Scintillator_box  = new G4Box("Scintillator_A", x4, y0,
                                      ScintillatorThickness);

  auto* Scintillator_pir = new G4GenericTrap("Scintillator_B", ScintillatorThickness,
                                             std::vector<G4TwoVector>{{0,-2*y0},
       {0,-2*y0}, {-x4,-y0}, {x4,-y0}, {0,-2*y0}, {0,-2*y0}, {-x4,-y0}, {x4,-y0}});

  auto* Scintillator_pir2 = new G4GenericTrap("Scintillator_C", ScintillatorThickness,
                                              std::vector<G4TwoVector>{{0,2*y0},
        {0,2*y0}, {x4,y0}, {-x4,y0}, {0,2*y0}, {0,2*y0}, {x4,y0}, {-x4,y0}});

  Scintillator_log_A
    = new G4LogicalVolume(Scintillator_box, ScintillatorMaterial, "Scint_box", 0, 0, 0);
  Scintillator_log_B     = new G4LogicalVolume(Scintillator_pir, ScintillatorMaterial, "Scint_pir", 0, 0, 0);
  Scintillator_log_C     = new G4LogicalVolume(Scintillator_pir2, ScintillatorMaterial, "Scint_pir2", 0, 0, 0);

    // ************* Placement
  Scintillator_log_A->SetVisAttributes(VisAtt_Scintillator_log);
  Scintillator_log_B->SetVisAttributes(VisAtt_Scintillator_log);
  Scintillator_log_C->SetVisAttributes(VisAtt_Scintillator_log);

  fScintillator_A = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), Scintillator_log_A,
                    "Scintillator_box", world_LV, false, 0);
  fScintillator_B = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                     Scintillator_log_B, "Scintillator_pir", world_LV, false, 0);
  fScintillator_C = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                      Scintillator_log_C, "Scintillator_pir2", world_LV, false, 0);
}

void DetectorConstruction::AddCScintillator(G4VisAttributes* VisAtt_Scintillator_log, G4LogicalVolume* world_LV,
                                              G4double ScintillatorLength, G4double ScintillatorThickness,
                                              G4int ScintillatorMaterial) {

   G4Tubs* Scintillator_tube
           = new G4Tubs ("Scintillator", 0*m, ScintillatorLength, ScintillatorThickness,
                         0*deg, 360*deg);
  if (ScintillatorMaterial == 0) {
    Scintillator_log
            = new G4LogicalVolume(Scintillator_tube, CeBr_3, "Scintillator", 0, 0, 0);
  }
  else if (ScintillatorMaterial ==1) {
    Scintillator_log
            = new G4LogicalVolume(Scintillator_tube, CsI_Tl, "Scintillator", 0, 0, 0);
  }
  else if (ScintillatorMaterial ==2) {
    Scintillator_log
            = new G4LogicalVolume(Scintillator_tube, NaI_Tl, "Scintillator", 0, 0, 0);
  }
 // ************* Placement
  Scintillator_log->SetVisAttributes(VisAtt_Scintillator_log);

  fScintillator
          = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Scintillator_log, "Scintillator",
                              world_LV, false, 0);
}

void DetectorConstruction::ChooseReflector(G4int ReflectorType, G4OpticalSurface* Sc_Ref, G4int ScintillatorMaterial) {
  Sc_Ref->SetType(dielectric_dielectric);
  Sc_Ref->SetModel(unified);
  G4double PhotonEnergy [nEntries];
  G4double rIndex [nEntries];

  if (ScintillatorMaterial == 0) { //CeBr3
    std::copy(std::begin(PhotonEnergyCeBr3), std::end(PhotonEnergyCeBr3), PhotonEnergy);
    std::copy(std::begin(rIndexCeBr3), std::end(rIndexCeBr3), rIndex);
  }
  else if (ScintillatorMaterial ==1) {  //CsI
    std::copy(std::begin(PhotonEnergyCsI), std::end(PhotonEnergyCsI), PhotonEnergy);
    std::copy(std::begin(rIndexCsI), std::end(rIndexCsI), rIndex);
  }
  else if (ScintillatorMaterial ==2) {  //CsI
    std::copy(std::begin(PhotonEnergyNaI), std::end(PhotonEnergyNaI), PhotonEnergy);
    std::copy(std::begin(rIndexNaI), std::end(rIndexNaI), rIndex);
  }

  if (ReflectorType == 0) {
    Sc_Ref->SetFinish(groundbackpainted);
    Sc_Ref->SetSigmaAlpha(0.0);
    fSc_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex, nEntries);
    fSc_RefMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity, nEntries);
    fSc_RefMPT->AddProperty("SPECULARSPIKECONSTANT", PhotonEnergy, SpSp, nEntries);
  }
  if (ReflectorType == 1) {
    Sc_Ref->SetFinish(groundbackpainted);
    Sc_Ref->SetSigmaAlpha(0.0);
    fSc_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex_air, nEntries);
    fSc_RefMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity, nEntries);
    fSc_RefMPT->AddProperty("SPECULARSPIKECONSTANT", PhotonEnergy, SpSp, nEntries);
  }
    if (ReflectorType == 2) {
    Sc_Ref->SetFinish(polishedbackpainted);
    Sc_Ref->SetSigmaAlpha(0.0);
    fSc_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex_air, nEntries);
    fSc_RefMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity, nEntries);
    fSc_RefMPT->AddProperty("SPECULARSPIKECONSTANT", PhotonEnergy, SpSp, nEntries);
    }

  if (ReflectorType == 3) {
    Sc_Ref->SetFinish(groundbackpainted);
    Sc_Ref->SetSigmaAlpha(0.6);
    fSc_RefMPT->AddProperty("RINDEX", PhotonEnergy, rIndex_air, nEntries);
    fSc_RefMPT->AddProperty("REFLECTIVITY", PhotonEnergy, reflectivity, nEntries);
    fSc_RefMPT->AddProperty("LAMBERTIAN", PhotonEnergy, SpSp, nEntries);
  }

  Sc_Ref->SetMaterialPropertiesTable(fSc_RefMPT);
}

