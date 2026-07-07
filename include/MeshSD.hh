#pragma once
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include <vector>
#include <string>

// ─────────────────────────────────────────────────────────────────
// MeshSD  – 3-D Cartesian mesh sensitive detector
//
// Replicates two FLUKA USRBIN cards on the SAME volume:
//
//  EdepSlab: ENERGY,  x[-100,100] y[-100,100] z[-THICK,0],  50×50×300 bins
//  NeutSlab: NEUTRON, x[-100,100] y[-100,100] z[-THICK,0], 100×100×105 bins
//
// A G4LogicalVolume can only carry ONE sensitive detector, so both grids
// (which use independent binnings) are scored by a single SD instance.
// ─────────────────────────────────────────────────────────────────

class MeshSD : public G4VSensitiveDetector
{
public:
    struct GridSpec {
        G4double xmin, xmax; G4int nx;
        G4double ymin, ymax; G4int ny;
        G4double zmin, zmax; G4int nz;
    };

    MeshSD(const G4String& name,
           const GridSpec& edepGrid,
           const GridSpec& neutGrid);
    ~MeshSD() override;

    void   Initialize(G4HCofThisEvent*) override {}
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;
    void   EndOfEvent(G4HCofThisEvent*) override {}

    void Reset();
    // Write Z-projections (XY-integrated) of both grids to CSV
    void Save(const std::string& tag, int nevents) const;

private:
    struct Grid {
        GridSpec spec;
        G4double dx, dy, dz;
        std::vector<double> data;
        int BinIndex(G4double x, G4double y, G4double z) const;
        void SaveZProjection(const std::string& fname, const std::string& particleTag,
                              double norm) const;
    };

    Grid fEdep;
    Grid fNeut;

    static Grid MakeGrid(const GridSpec& spec);
};
