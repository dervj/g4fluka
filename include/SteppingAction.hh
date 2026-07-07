#pragma once
#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>
#include <string>

class DetectorConstruction;
class G4Step;

// ─────────────────────────────────────────────────────────────────
// SteppingAction
//
// Scores particles crossing the SLAB→DSTRM boundary (z = 0 plane),
// replicating the FLUKA USRBDX cards:
//   BdxNeut  – neutrons,  0 → 500 MeV, 100 log-E bins
//   BdxProt  – protons,   0.1 keV → 500 MeV, 100 log-E bins
//   BdxPhot  – photons,   0.1 keV → 500 MeV, 100 log-E bins
//   BdxCHad  – charged hadrons (pions±, kaons±, etc.)
//   BdxElec  – e+/e-
//
// Also scores the double-differential neutron yield (USRYIELD):
//   YldNeut: d²N/dEdΩ, 18 polar-angle bins [0°,180°], kinetic energy
//            from 1e-9 GeV to 0.5 GeV.
//
// All data are accumulated per run and written to ROOT TH1D/TH2D,
// and also to plain-text CSV files for easy comparison with FLUKA.
// ─────────────────────────────────────────────────────────────────

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(const DetectorConstruction* det);
    ~SteppingAction() override;

    void UserSteppingAction(const G4Step* step) override;

    // Called by RunAction
    void Reset();
    void Save(const std::string& tag);   // tag = run label, e.g. "W_500MeV"

private:
    const DetectorConstruction* fDet;

    // Energy bin edges (log-spaced, built once)
    std::vector<double> fNeutBins;   // 101 edges, 1e-9 GeV → 0.5 GeV
    std::vector<double> fCharBins;   // 101 edges, 1e-4 GeV → 0.5 GeV
    std::vector<double> fAnglBins;   // 19 edges,  0° → 180°

    // Accumulators [iAngle][iEnergy]  – for YldNeut
    std::vector<std::vector<double>> fYldNeut;
    // 1-D spectra [iEnergy]
    std::vector<double> fBdxNeut;
    std::vector<double> fBdxProt;
    std::vector<double> fBdxPhot;
    std::vector<double> fBdxCHad;
    std::vector<double> fBdxElec;

    // Helpers
    int  FindBin(const std::vector<double>& edges, double val) const;
    void WriteCSV(const std::string& filename,
                  const std::vector<double>& edges,
                  const std::vector<double>& counts,
                  double norm) const;
    void Write2DCSV(const std::string& filename,
                    const std::vector<double>& eEdges,
                    const std::vector<double>& aEdges,
                    const std::vector<std::vector<double>>& counts,
                    double norm) const;

    G4double fSlabExitZ = 0.;   // z position of slab downstream face (cm)
    G4int    fNevents   = 0;
};
