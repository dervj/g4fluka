#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

// ─────────────────────────────────────────────────────────────────
// DetectorConstruction
//
// Geometry mirrors the FLUKA g4fluka3.inp exactly:
//
//   World (vacuum)     : x[-100,100], y[-100,100], z[-55,200]  cm
//   UPSTRM (vacuum)    : same XY,     z[-55, 0]   cm  (minus slab)
//   SLAB   (material)  : same XY,     z[-THICK, 0] cm
//   DSTRM  (vacuum)    : same XY,     z[0, 200]   cm  (minus slab)
//
// Beam travels along +Z, source at z = -53 cm.
// THICK = 50 cm  → Tungsten (default)
// THICK = 10.5 cm → BPE  (set via /det/setMaterial and /det/setThickness)
// ─────────────────────────────────────────────────────────────────

class G4Material;
class G4VPhysicalVolume;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Called by messenger or main()
    void SetSlabMaterial(const G4String& name);
    void SetSlabThickness(G4double thick);   // in cm (will be converted internally)

    G4double       GetSlabThickness()  const { return fSlabThick; }
    G4Material*    GetSlabMaterial()   const { return fSlabMat;   }
    G4LogicalVolume* GetSlabLV()       const { return fSlabLV;    }
    G4LogicalVolume* GetDstrmLV()      const { return fDstrmLV;   }

private:
    void DefineMaterials();

    G4double  fSlabThick;     // half-thickness for G4Box (cm stored, converted)
    G4Material* fSlabMat;
    G4Material* fMatTungsten;
    G4Material* fMatBPE;
    G4Material* fMatVacuum;

    G4LogicalVolume* fSlabLV  = nullptr;
    G4LogicalVolume* fDstrmLV = nullptr;

    DetectorMessenger* fMessenger = nullptr;
};
