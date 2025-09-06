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
/// \file optical/OpNovice2/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4NistManager.hh"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"



class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(
        G4int detectorForm=1, // 0-cyl, 1-quad, 2-hexagonal
        G4bool hasReflector=true,
        G4bool hasAlu=true,
        G4double ScintillatorThickness=10 *mm,
        G4int ReflectorType=0,
        G4double ScintillatorSideLength=51 *mm,
        G4double ScintillatorRadius=25.5 *mm,
        G4int ScintillatorMaterial=0, //CeBr3
        G4double ScintillatorSideWidth=51 *mm,
        G4int NumberOfSiPM=1
    );

    virtual ~DetectorConstruction();

    G4VPhysicalVolume* GetScintillator() {return fScintillator;}
    G4double GetTotalThickness() const {
        if (fHasReflector) {
          if (fHasAlu) { return fScintillatorThickness + 2*fSiPM_z; }
          else { return fScintillatorThickness + 2*fRefl_z; }
        }
        else {return fScintillatorThickness;}
      }


    //G4OpticalSurface* GetSurface(void) {return fSurface;}

    void SetSurfaceFinish(const G4OpticalSurfaceFinish finish) {
      fSurface->SetFinish(finish);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    G4OpticalSurfaceFinish GetSurfaceFinish(void)
      {return fSurface->GetFinish();}

    void SetSurfaceType(const G4SurfaceType type) {
      fSurface->SetType(type);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }

    void SetSurfaceModel(const G4OpticalSurfaceModel model) {
      fSurface->SetModel(model);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    G4OpticalSurfaceModel GetSurfaceModel(void)
      {return fSurface->GetModel();}

    void SetSurfaceSigmaAlpha(G4double v);

    void SetDetectorSection(G4int v);

    void SetDetectorHasReflector(G4bool v);

    void SetDetectorHasAlu(G4bool v);

    void ChangeNumberOfSiPM(G4int v);

    void SetCrystalThickness(G4double v);

    void SetReflectorType(G4int v);

    void SetQuadraticCrystalSideLength(G4double v);
    void SetQuadraticCrystalWidth(G4double v);

    void SetCylindricalCrystalRadius(G4double v);

    void SetScintillatorMaterial(G4int v);

    void AddBoxMPV(const char* c, G4MaterialPropertyVector* mpv);
    void AddBoxMPCV(const char* c, G4double v);
    G4MaterialPropertiesTable* GetBoxMaterialPropertiesTable()
      {return fBoxMPT;}

    void AddWorldMPV(const char* c, G4MaterialPropertyVector* mpv);
    void AddWorldMPCV(const char* c, G4double v);
    G4MaterialPropertiesTable* GetWorldMaterialPropertiesTable()
      {return fWorldMPT;}

    void AddSurfaceMPV(const char* c, G4MaterialPropertyVector* mpv);
    G4MaterialPropertiesTable* GetSurfaceMaterialPropertiesTable()
      {return fSurfaceMPT;}



    virtual G4VPhysicalVolume* Construct();

  private:
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

//    G4VPhysicalVolume* fTank;     //G4double fTank_x;     //G4double fTank_y; // G4double fTank_z;

    G4VPhysicalVolume* fScintillator;
    G4VPhysicalVolume* fScintillator_A;
    G4VPhysicalVolume* fScintillator_B;
    G4VPhysicalVolume* fScintillator_C;

    G4double fScintillatorThickness;
    G4double fScintillatorLength;
    G4double fScintillatorWidth;
    G4double fScintillatorRadius;
    G4int fScintillatorMaterial;

    G4LogicalVolume* Scintillator_log;
    G4LogicalVolume* Scintillator_log_A;
    G4LogicalVolume* Scintillator_log_B;
    G4LogicalVolume* Scintillator_log_C;

    G4LogicalVolume* AluBoxLogXY;
    G4LogicalVolume* AluBoxLogXZ;
    G4LogicalVolume* AluBoxLogYZ;
    G4LogicalVolume* AluBoxLog;
    G4LogicalVolume* AluTrapLog1, * AluTrapLog2, * AluSideTrapLog1, * AluSideTrapLog2;


    G4LogicalVolume* SiPMLogXY;
    G4LogicalVolume* SiPMLogYZ;
    G4LogicalVolume* SiPMLogTrap01;
    G4LogicalVolume* SiPMLogTrap02;
    G4LogicalVolume* SiPMLogTrap11;
    G4LogicalVolume* SiPMLogTrap12;


    G4int fNumberOfSiPM;
    G4int fDetectorSection;
    G4bool fHasReflector;
    G4bool fHasAlu;
    G4int  fReflectorType;

    G4VPhysicalVolume* fSiPMXY;
    G4VPhysicalVolume* fSiPMYZ2;
    G4VPhysicalVolume* fSiPMYZ1;

    G4VPhysicalVolume* fSiPMTrap01;
    G4VPhysicalVolume* fSiPMTrap02;
    G4VPhysicalVolume* fSiPMTrap11;
    G4VPhysicalVolume* fSiPMTrap12;

    G4VPhysicalVolume* fReflectorYZ1;
    G4VPhysicalVolume* fReflectorYZ2;
    G4VPhysicalVolume* fReflectorXZ1;
    G4VPhysicalVolume* fReflectorXZ2;
    G4VPhysicalVolume* fReflectorXY1;
    G4VPhysicalVolume* fReflector;
    G4VPhysicalVolume* fReflector_2;
    G4VPhysicalVolume* fAlu;
    G4VPhysicalVolume* fAlu_3;
    G4VPhysicalVolume* fReflectorBox1;
    G4VPhysicalVolume* fReflectorBox2;
    G4VPhysicalVolume* fReflectorS1;
    G4VPhysicalVolume* fReflectorS2;
    G4VPhysicalVolume* fReflectorTrap21;
    G4VPhysicalVolume* fReflectorTrap22;
    G4VPhysicalVolume* fReflectorTrap11;
    G4VPhysicalVolume* fReflectorTrap12;

    G4VPhysicalVolume* fAluBoxXY1;
    G4VPhysicalVolume* fAluBoxXZ1;
    G4VPhysicalVolume* fAluBoxXZ2;
    G4VPhysicalVolume* fAluBoxYZ1;
    G4VPhysicalVolume* fAluBoxYZ2;

    G4VPhysicalVolume* fAluBox1;
    G4VPhysicalVolume* fAluBox2;
    G4VPhysicalVolume* fAluTrap11, * fAluTrap12, * fAluTrap21, * fAluTrap22, * fAluSideTrap1, * fAluSideTrap2;


    G4LogicalBorderSurface* Sc_RefLogYZ1; //quadratic
    G4LogicalBorderSurface* Sc_RefLogYZ2;  //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_21; //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_22;  //quadratic
    G4LogicalBorderSurface* Sc_Ref_log_3;  //quadratic

    G4LogicalBorderSurface* Sc_Ref_log; //cylindrical
    G4LogicalBorderSurface* Sc_Ref_log_2; //cylindrical

    G4LogicalBorderSurface* Sc_RefLogTr11; //H
    G4LogicalBorderSurface* Sc_RefLogTr12;  //H
    G4LogicalBorderSurface* Sc_RefLogTr21; //H
    G4LogicalBorderSurface* Sc_RefLogTr22;  //H
    G4LogicalBorderSurface* Sc_RefLogBox1;  //H
    G4LogicalBorderSurface* Sc_RefLogBox2;  //H
    G4LogicalBorderSurface* Sc_RefLogSide1;  //H
    G4LogicalBorderSurface* Sc_RefLogSide2;  //H


    G4double fSiPM_r;
    G4double fSiPM_z;
    G4double fRefl_z;
    G4double fAlu_thickness;

    G4VPhysicalVolume* fReflector2;

   // G4Material* Vacuum;

    G4double fRef_r1;
    G4double fRef_r2;

    G4double z;
    G4double a;

    G4double y0, y1, y2, y3, y4, y6;
    G4double x0, x1, x2, x3, x4, x5;


    G4OpticalSurface* fSurface;

    DetectorMessenger* fDetectorMessenger;

    G4MaterialPropertiesTable* fBoxMPT;
    G4MaterialPropertiesTable* fReflectorMPT;
    G4MaterialPropertiesTable* fWorldMPT;
    G4MaterialPropertiesTable* fScintillatorMPT;
    G4MaterialPropertiesTable* fSurfaceMPT;
    G4MaterialPropertiesTable* fSiPMMPT;
    G4MaterialPropertiesTable* fSc_RefMPT;
    G4MaterialPropertiesTable* fSc_SiPMMPT;
    G4MaterialPropertiesTable* fSi_RefMPT;

    void DefineMaterials();
    void DefineOpticalProperties(G4int ScintillatorMaterial);
    void AddQHousing(G4VisAttributes* VisAtt_Alu_log, G4LogicalVolume* world_LV, G4int NumberOfSiPM, G4double AluLength, G4double AluWidth, G4double Thickness);
    void AddHHousing(G4VisAttributes* VisAtt_Alu_log, G4LogicalVolume* world_LV, G4double ScintillatorLength,
                     G4double Thickness, G4double ScintillatorThickness);

    void AddQSiPM(G4VisAttributes* VisAtt_SiPM_log, G4LogicalVolume* world_LV, G4int NumberOfSiPM, G4double ScintXSize,
                 G4double ScintYSize, G4double ScintZSize, G4double SiPMThickness);
    void AddHSiPM(G4VisAttributes* VisAtt_SiPM_log, G4LogicalVolume* world_LV, G4double ScintillatorLength,
                  G4double ScintillatorThickness, G4double SiPMThickness);

    void AddCSiPM(G4VisAttributes* VisAtt_SiPM_log, G4LogicalVolume* world_LV,
                 G4double SiPMSize, G4double SiPM_z);

    void AddQReflector(G4VisAttributes* VisAtt_Reflector_log, G4LogicalVolume* world_LV, G4int NumberOfSiPM,
                       G4double ScintillatorSideLength, G4double ScintillatorSideWidth, G4double SiPM_z, G4double Refl_z);
    void AddHReflector(G4VisAttributes* VisAtt_Reflector_log, G4LogicalVolume* world_LV,
                       G4double ScintillatorSideLength, G4double ScintillatorThickness, G4double Refl_z);

    void AddCHousing(G4VisAttributes* VisAtt_Alu_log, G4LogicalVolume* world_LV, G4double AluLength, G4double SiPM_z);
    void AddCReflector(G4VisAttributes* VisAtt_Reflector_log, G4LogicalVolume* world_LV,
                       G4double ScintillatorRadius, G4double SiPM_z, G4double Refl_z);
    void AddQScintillator(G4VisAttributes* VisAtt_Scintillator_log, G4LogicalVolume* world_LV,
                         G4double ScintillatorLength, G4double ScintillatorWidth, G4double ScintillatorThickness, G4int ScintillatorMaterial);
    void AddHScintillator(G4VisAttributes* VisAtt_Scintillator_log, G4LogicalVolume* world_LV,
                          G4double ScintillatorLength, G4double ScintillatorThickness, G4Material* ScintillatorMaterial);

    void AddCScintillator(G4VisAttributes* VisAtt_Scintillator_log, G4LogicalVolume* world_LV,
                         G4double ScintillatorLength, G4double ScintillatorThickness, G4int ScintillatorMaterial);

    void ChooseReflector(G4int ReflectorType, G4OpticalSurface* Sc_Ref, G4int ScintillatorMaterial);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
