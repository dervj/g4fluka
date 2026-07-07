#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sstream>

// ─── Helper: log-spaced bin edges ────────────────────────────────
static std::vector<double> LogEdges(double emin_GeV, double emax_GeV, int nbins)
{
    std::vector<double> v(nbins+1);
    double step = std::log10(emax_GeV/emin_GeV)/nbins;
    for(int i=0;i<=nbins;++i)
        v[i] = emin_GeV * std::pow(10., i*step);
    return v;
}

// ─── Helper: linear bin edges ────────────────────────────────────
static std::vector<double> LinEdges(double lo, double hi, int nbins)
{
    std::vector<double> v(nbins+1);
    double step = (hi-lo)/nbins;
    for(int i=0;i<=nbins;++i) v[i]=lo+i*step;
    return v;
}

// ─────────────────────────────────────────────────────────────────
SteppingAction::SteppingAction(const DetectorConstruction* det)
    : fDet(det)
{
    // FLUKA USRBDX BdxNeut: Emin=1e-9 GeV, Emax=0.5 GeV, 100 log-bins
    fNeutBins = LogEdges(1.e-9, 0.5, 100);
    // FLUKA USRBDX BdxProt/Phot/CHad/Elec: Emin=1e-4 GeV, Emax=0.5 GeV, 100 log-bins
    fCharBins = LogEdges(1.e-4, 0.5, 100);
    // FLUKA USRYIELD YldNeut: 18 angle bins [0,180] degrees, linear
    fAnglBins = LinEdges(0., 180., 18);

    Reset();
}

SteppingAction::~SteppingAction() {}

void SteppingAction::Reset()
{
    fBdxNeut.assign(100, 0.);
    fBdxProt.assign(100, 0.);
    fBdxPhot.assign(100, 0.);
    fBdxCHad.assign(100, 0.);
    fBdxElec.assign(100, 0.);

    fYldNeut.assign(18, std::vector<double>(100, 0.));
    // (for yield: 18 angle bins × 100 energy bins matching BdxNeut range)

    fNevents = 0;
}

// ─────────────────────────────────────────────────────────────────
int SteppingAction::FindBin(const std::vector<double>& edges, double val) const
{
    if (val < edges.front() || val >= edges.back()) return -1;
    auto it = std::upper_bound(edges.begin(), edges.end(), val);
    return (int)(it - edges.begin()) - 1;
}

// ─────────────────────────────────────────────────────────────────
void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // Only score at boundary crossings
    G4StepPoint* postPt = step->GetPostStepPoint();
    if (postPt->GetStepStatus() != fGeomBoundary) return;

    G4StepPoint* prePt = step->GetPreStepPoint();
    const G4VPhysicalVolume* preVol  = prePt->GetPhysicalVolume();
    const G4VPhysicalVolume* postVol = postPt->GetPhysicalVolume();
    if (!preVol || !postVol) return;

    // We want crossings from SLABREG (name "SLABREG") → DSTREG (name "DSTREG")
    bool fromSlab = (preVol->GetName()  == "SLABREG");
    bool toDstrm  = (postVol->GetName() == "DSTREG");
    if (!fromSlab || !toDstrm) return;

    const G4Track* track = step->GetTrack();
    const G4ParticleDefinition* pdef = track->GetParticleDefinition();
    G4double ekin = track->GetKineticEnergy();
    G4double eGeV = ekin / GeV;

    G4ThreeVector dir = track->GetMomentumDirection();
    // Polar angle w.r.t. +Z (beam direction)
    G4double theta_deg = std::acos(std::clamp(dir.z(), -1.0, 1.0)) * (180./M_PI);

    // ── BdxNeut ──────────────────────────────────────────────────
    if (pdef == G4Neutron::Definition()) {
        int ib = FindBin(fNeutBins, eGeV);
        if (ib >= 0) {
            fBdxNeut[ib] += 1.;
            // YldNeut: d²N/dEdΩ  (same energy binning as BdxNeut)
            int ia = FindBin(fAnglBins, theta_deg);
            if (ia >= 0) fYldNeut[ia][ib] += 1.;
        }
    }
    // ── BdxProt ──────────────────────────────────────────────────
    else if (pdef == G4Proton::Definition()) {
        int ib = FindBin(fCharBins, eGeV);
        if (ib >= 0) fBdxProt[ib] += 1.;
    }
    // ── BdxPhot ──────────────────────────────────────────────────
    else if (pdef == G4Gamma::Definition()) {
        int ib = FindBin(fCharBins, eGeV);
        if (ib >= 0) fBdxPhot[ib] += 1.;
    }
    // ── BdxCHad (charged hadrons: π±, K±, etc.) ─────────────────
    else if (pdef->GetPDGCharge() != 0 &&
             pdef->GetParticleType() == "baryon" &&
             pdef != G4Proton::Definition()) {
        int ib = FindBin(fCharBins, eGeV);
        if (ib >= 0) fBdxCHad[ib] += 1.;
    }
    else if (pdef->GetPDGCharge() != 0 &&
             pdef->GetParticleType() == "meson") {
        int ib = FindBin(fCharBins, eGeV);
        if (ib >= 0) fBdxCHad[ib] += 1.;
    }
    // ── BdxElec (e+/e-) ──────────────────────────────────────────
    else if (pdef == G4Electron::Definition() ||
             pdef == G4Positron::Definition()) {
        int ib = FindBin(fCharBins, eGeV);
        if (ib >= 0) fBdxElec[ib] += 1.;
    }
}

// ─────────────────────────────────────────────────────────────────
// Save results: normalise per primary, per cm² (area=40000 cm²),
//               per GeV (bin width), consistent with FLUKA USRBDX Φ1
// ─────────────────────────────────────────────────────────────────
void SteppingAction::Save(const std::string& tag)
{
    // Get number of events from run manager
    const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
    double N = run ? (double)run->GetNumberOfEvent() : 1.;
    double area_cm2 = 40000.;   // 200×200 cm² (FLUKA USRBDX Area)

    WriteCSV("BdxNeut_" + tag + ".csv", fNeutBins, fBdxNeut, N * area_cm2);
    WriteCSV("BdxProt_" + tag + ".csv", fCharBins, fBdxProt, N * area_cm2);
    WriteCSV("BdxPhot_" + tag + ".csv", fCharBins, fBdxPhot, N * area_cm2);
    WriteCSV("BdxCHad_" + tag + ".csv", fCharBins, fBdxCHad, N * area_cm2);
    WriteCSV("BdxElec_" + tag + ".csv", fCharBins, fBdxElec, N * area_cm2);
    Write2DCSV("YldNeut_" + tag + ".csv", fNeutBins, fAnglBins, fYldNeut, N);
}

void SteppingAction::WriteCSV(const std::string& fname,
                               const std::vector<double>& edges,
                               const std::vector<double>& counts,
                               double norm) const
{
    std::ofstream f(fname);
    f << "# Elow_GeV,Ehigh_GeV,dPhi_dE_per_primary_per_cm2_per_GeV\n";
    for (int i=0; i<(int)counts.size(); ++i) {
        double dE = edges[i+1] - edges[i];
        double val = (norm > 0 && dE > 0) ? counts[i]/(norm*dE) : 0.;
        f << std::scientific << std::setprecision(6)
          << edges[i] << "," << edges[i+1] << "," << val << "\n";
    }
    f.close();
    G4cout << "[SteppingAction] Wrote " << fname << G4endl;
}

void SteppingAction::Write2DCSV(const std::string& fname,
                                 const std::vector<double>& eEdges,
                                 const std::vector<double>& aEdges,
                                 const std::vector<std::vector<double>>& counts,
                                 double norm) const
{
    std::ofstream f(fname);
    f << "# theta_lo_deg,theta_hi_deg,Elow_GeV,Ehigh_GeV,d2N_dEdOmega_per_primary\n";
    for (int ia=0; ia<(int)counts.size(); ++ia) {
        // solid angle of ring: dΩ = 2π sin(θ) dθ
        double tlo = aEdges[ia]   * M_PI/180.;
        double thi = aEdges[ia+1] * M_PI/180.;
        double dOmega = 2.*M_PI*(std::cos(tlo) - std::cos(thi));
        for (int ie=0; ie<(int)counts[ia].size(); ++ie) {
            double dE = eEdges[ie+1] - eEdges[ie];
            double val = (norm > 0 && dE > 0 && std::abs(dOmega) > 0)
                         ? counts[ia][ie]/(norm*dE*std::abs(dOmega))
                         : 0.;
            f << std::scientific << std::setprecision(6)
              << aEdges[ia] << "," << aEdges[ia+1] << ","
              << eEdges[ie] << "," << eEdges[ie+1] << ","
              << val << "\n";
        }
    }
    f.close();
    G4cout << "[SteppingAction] Wrote " << fname << G4endl;
}
