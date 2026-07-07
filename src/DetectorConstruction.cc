#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "MeshSD.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"

#include <algorithm>
#include <cmath>

// ─── Constructor ──────────────────────────────────────────────────
DetectorConstruction::DetectorConstruction()
    : fSlabThick(50.*cm),   // default: 50 cm Tungsten  (THICK=50 in .inp)
      fSlabMat(nullptr), fMatTungsten(nullptr), fMatBPE(nullptr), fMatVacuum(nullptr)
{
    fMessenger = new DetectorMessenger(this);
    DefineMaterials();
    fSlabMat = fMatTungsten;
}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

// ─── Materials ────────────────────────────────────────────────────
// Matches FLUKA:
//   BORON        Z=5,  rho=2.37 g/cm3  (elemental)
//   POLYETHY     rho=0.94, C2H4 by mass fraction
//   BPE          rho=0.95, 5% Boron + 85.6% Polyethylene + 9.4% Hydrogen
//   TUNGSTEN     pre-defined (NIST)
void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    // ── Vacuum ──────────────────────────────────────────────────
    fMatVacuum = nist->FindOrBuildMaterial("G4_Galactic");

    // ── Tungsten ────────────────────────────────────────────────
    fMatTungsten = nist->FindOrBuildMaterial("G4_W");

    // ── Elemental Boron  (FLUKA: MATERIAL Z=5, rho=2.37) ────────
    G4Material* matBoron = new G4Material("FlukaBoron",
                                          /*Z=*/5., /*A=*/10.811*g/mole,
                                          2.37*g/cm3);

    // ── Polyethylene  (C2H4)n  rho=0.94  ────────────────────────
    // FLUKA COMPOUND: -0.143711 HYDROGEN -0.856289 CARBON
    G4Element* elH = nist->FindOrBuildElement("H");
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Material* matPoly = new G4Material("FlukaPolyethy",
                                          0.94*g/cm3, /*nComponents=*/2);
    matPoly->AddElement(elH, 0.143711);
    matPoly->AddElement(elC, 0.856289);

    // ── BPE  rho=0.95  ───────────────────────────────────────────
    // FLUKA COMPOUND: -0.05 BORON  -0.856 POLYETHY  -0.094 HYDROGEN
    // Note: FLUKA negative fractions = mass fractions
    fMatBPE = new G4Material("FlukaBPE", 0.95*g/cm3, /*nComponents=*/3);
    fMatBPE->AddMaterial(matBoron,  0.050);
    fMatBPE->AddMaterial(matPoly,   0.856);
    fMatBPE->AddElement (elH,       0.094);
}

// ─── Geometry ────────────────────────────────────────────────────
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // ─── World box: x[-100,100] y[-100,100], covering z[-55,200] cm ───
    // The world volume must be placed at the origin with no offset (its
    // frame IS the global frame), so it has to be symmetric about z=0.
    // Daughters below use their true absolute-frame z coordinates.
    G4double wx = 100.*cm, wy = 100.*cm;
    G4double wzMin = -55.*cm, wzMax = 200.*cm;
    G4double wzHalf = std::max(std::abs(wzMin), std::abs(wzMax));

    G4Box* worldBox = new G4Box("World", wx, wy, wzHalf);
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldBox, fMatVacuum, "World");
    G4VPhysicalVolume* worldPV = new G4PVPlacement(
        nullptr, G4ThreeVector(0,0,0), worldLV, "World", nullptr, false, 0);

    // ─── Slab: x[-100,100] y[-100,100] z[-THICK, 0] ─────────────
    G4double slabHalf = fSlabThick/2.;
    G4double slabZ    = -slabHalf;   // centre in world coords

    G4Box* slabBox = new G4Box("Slab", wx, wy, slabHalf);
    fSlabLV = new G4LogicalVolume(slabBox, fSlabMat, "SlabLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,slabZ), fSlabLV,
                      "SLABREG", worldLV, false, 0, true);

    // ─── Downstream volume: z[0, 200] cm ─────────────────────────
    G4double dHalf   = 100.*cm;
    G4double dCentre = 100.*cm;   // [0,200] in world frame
    G4Box* dBox = new G4Box("Dstrm", wx, wy, dHalf);
    fDstrmLV = new G4LogicalVolume(dBox, fMatVacuum, "DstrmLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,dCentre), fDstrmLV,
                      "DSTREG", worldLV, false, 0, true);

    // ─── Upstream volume fills the rest (world - slab - dstrm) ───
    // Represented implicitly by the world volume.
    // (No explicit placement needed; Geant4 tracks in the world.)

    // ── Visualisation ────────────────────────────────────────────
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4VisAttributes* slabVis = new G4VisAttributes(G4Colour(0.8,0.8,0.0,0.4));
    slabVis->SetForceSolid(true);
    fSlabLV->SetVisAttributes(slabVis);
    G4VisAttributes* dVis = new G4VisAttributes(G4Colour(0.4,0.7,1.0,0.2));
    dVis->SetForceSolid(false);
    fDstrmLV->SetVisAttributes(dVis);

    return worldPV;
}

// ─── Sensitive Detectors ─────────────────────────────────────────
void DetectorConstruction::ConstructSDandField()
{
    G4SDManager* sdMgr = G4SDManager::GetSDMpointer();

    // A G4LogicalVolume can only carry one SD, so both USRBIN-equivalent
    // meshes (EdepSlab: 50x50x300, NeutSlab: 100x100x105) are scored by a
    // single MeshSD instance attached once to the slab.
    MeshSD::GridSpec edepSpec{-100.*cm, 100.*cm, 50,
                               -100.*cm, 100.*cm, 50,
                               -fSlabThick, 0.*cm, 300};
    MeshSD::GridSpec neutSpec{-100.*cm, 100.*cm, 100,
                               -100.*cm, 100.*cm, 100,
                               -fSlabThick, 0.*cm, 105};
    auto* sdSlab = new MeshSD("SlabMesh_SD", edepSpec, neutSpec);
    sdMgr->AddNewDetector(sdSlab);
    fSlabLV->SetSensitiveDetector(sdSlab);
}

// ─── Setters called by messenger ─────────────────────────────────
void DetectorConstruction::SetSlabMaterial(const G4String& name)
{
    if (name == "BPE" || name == "bpe") {
        fSlabMat = fMatBPE;
    } else {
        fSlabMat = fMatTungsten;
    }
    if (fSlabLV) {
        fSlabLV->SetMaterial(fSlabMat);
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
}

void DetectorConstruction::SetSlabThickness(G4double thick)
{
    fSlabThick = thick;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}
