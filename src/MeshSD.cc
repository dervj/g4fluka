#include "MeshSD.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Neutron.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iomanip>
#include <numeric>
#include <cmath>

MeshSD::Grid MeshSD::MakeGrid(const GridSpec& spec)
{
    Grid g;
    g.spec = spec;
    g.dx = (spec.xmax - spec.xmin)/spec.nx;
    g.dy = (spec.ymax - spec.ymin)/spec.ny;
    g.dz = (spec.zmax - spec.zmin)/spec.nz;
    g.data.assign(spec.nx*spec.ny*spec.nz, 0.);
    return g;
}

void MeshSD::Grid::AddAlongZ(G4double x, G4double y, G4double z1, G4double z2, G4double totalValue)
{
    int ix = (int)std::floor((x - spec.xmin)/dx);
    int iy = (int)std::floor((y - spec.ymin)/dy);
    if (ix<0||ix>=spec.nx||iy<0||iy>=spec.ny) return;

    G4double zLo = std::min(z1, z2);
    G4double zHi = std::max(z1, z2);
    G4double dzTotal = zHi - zLo;

    // Zero (or negligible) z displacement: can't apportion by overlap,
    // so credit it all to the single bin containing the step.
    if (dzTotal <= 0.) {
        int iz = (int)std::floor((zLo - spec.zmin)/dz);
        if (iz>=0 && iz<spec.nz)
            data[iz*spec.ny*spec.nx + iy*spec.nx + ix] += totalValue;
        return;
    }

    int izLo = std::max(0, (int)std::floor((zLo - spec.zmin)/dz));
    int izHi = std::min(spec.nz-1, (int)std::floor((zHi - spec.zmin)/dz));

    for (int iz = izLo; iz <= izHi; ++iz) {
        G4double binLo = spec.zmin + iz*dz;
        G4double binHi = binLo + dz;
        G4double overlap = std::min(zHi, binHi) - std::max(zLo, binLo);
        if (overlap <= 0.) continue;
        data[iz*spec.ny*spec.nx + iy*spec.nx + ix] += totalValue * (overlap/dzTotal);
    }
}

void MeshSD::Grid::SaveZProjection(const std::string& fname, const std::string& /*particleTag*/,
                                    double norm) const
{
    std::vector<double> zProf(spec.nz, 0.);
    for (int iz=0; iz<spec.nz; ++iz) {
        for (int iy=0; iy<spec.ny; ++iy)
            for (int ix=0; ix<spec.nx; ++ix)
                zProf[iz] += data[iz*spec.ny*spec.nx + iy*spec.nx + ix];
    }

    std::ofstream f(fname);
    f << "# z_lo_cm, z_hi_cm, value_per_primary\n";
    for (int iz=0; iz<spec.nz; ++iz) {
        double zlo = spec.zmin + iz*dz;
        double val = zProf[iz] / norm;
        f << std::scientific << std::setprecision(6)
          << zlo/cm << "," << (zlo+dz)/cm << "," << val << "\n";
    }
    f.close();
    G4cout << "[MeshSD] Wrote " << fname << G4endl;
}

MeshSD::MeshSD(const G4String& name,
               const GridSpec& edepGrid,
               const GridSpec& neutGrid)
    : G4VSensitiveDetector(name),
      fEdep(MakeGrid(edepGrid)),
      fNeut(MakeGrid(neutGrid))
{
}

MeshSD::~MeshSD() {}

void MeshSD::Reset()
{
    std::fill(fEdep.data.begin(), fEdep.data.end(), 0.);
    std::fill(fNeut.data.begin(), fNeut.data.end(), 0.);
}

G4bool MeshSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4Track* track = step->GetTrack();
    G4ThreeVector prePos  = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();

    // Energy deposition in GeV (FLUKA USRBIN ENERGY), spread along the step's
    // z-extent rather than dumped entirely in the pre-step voxel.
    G4double edep = step->GetTotalEnergyDeposit() / GeV;
    if (edep > 0.) {
        fEdep.AddAlongZ(prePos.x(), prePos.y(), prePos.z(), postPos.z(), edep);
    }

    // Neutron track-length flux: tl/V  (FLUKA USRBIN NEUTRON), likewise
    // apportioned across every z bin the step actually crosses.
    if (track->GetParticleDefinition() == G4Neutron::Definition()) {
        G4double stepLen = step->GetStepLength() / cm;  // cm
        if (stepLen > 0.) {
            G4double vol = (fNeut.dx * fNeut.dy * fNeut.dz) / (cm*cm*cm);  // cm^3
            fNeut.AddAlongZ(prePos.x(), prePos.y(), prePos.z(), postPos.z(), stepLen / vol);
        }
    }

    return true;
}

void MeshSD::Save(const std::string& tag, int nevents) const
{
    double N = (nevents > 0) ? nevents : 1.;
    fEdep.SaveZProjection("Edep_Zproj_" + tag + ".csv", "Edep", N);
    fNeut.SaveZProjection("NeutFlux_Zproj_" + tag + ".csv", "NeutFlux", N);
}
