#include "RunAction.hh"
#include "SteppingAction.hh"
#include "MeshSD.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UserSteppingAction.hh"
#include "G4Material.hh"

#include <sstream>

// We access the SteppingAction via the RunManager's user action
static SteppingAction* GetSteppingAction()
{
    return dynamic_cast<SteppingAction*>(
        const_cast<G4UserSteppingAction*>(
            G4RunManager::GetRunManager()->GetUserSteppingAction()));
}

// Build a tag from the run's actual configuration (material, thickness,
// beam energy) rather than relying on an external env var — Geant4 has
// no built-in "/control/setenv" macro command to set one from a .mac file.
static std::string BuildRunTag()
{
    auto* rm = G4RunManager::GetRunManager();
    auto* det = dynamic_cast<const DetectorConstruction*>(rm->GetUserDetectorConstruction());
    auto* gen = dynamic_cast<const PrimaryGeneratorAction*>(rm->GetUserPrimaryGeneratorAction());

    std::string matName = (det && det->GetSlabMaterial()) ? det->GetSlabMaterial()->GetName() : "Mat";
    G4double thickCm = det ? det->GetSlabThickness()/CLHEP::cm : 0.;
    G4double energyMeV = gen ? gen->GetBeamEnergy()/CLHEP::MeV : 0.;

    std::ostringstream tag;
    tag << matName << "_" << thickCm << "cm_" << energyMeV << "MeV";
    return tag.str();
}

RunAction::RunAction() {}
RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
    SteppingAction* sa = GetSteppingAction();
    if (sa) sa->Reset();

    // Scan macros issue multiple /run/beamOn calls in one process; the mesh
    // must be cleared each time or successive energy points accumulate.
    auto* sdSlab = dynamic_cast<MeshSD*>(
        G4SDManager::GetSDMpointer()->FindSensitiveDetector("SlabMesh_SD", false));
    if (sdSlab) sdSlab->Reset();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    std::string tag = BuildRunTag();

    SteppingAction* sa = GetSteppingAction();
    if (sa) sa->Save(tag);

    // Save mesh SD (scores both Edep and NeutFlux grids)
    G4SDManager* sdMgr = G4SDManager::GetSDMpointer();
    int nev = run->GetNumberOfEvent();

    auto* sdSlab = dynamic_cast<MeshSD*>(sdMgr->FindSensitiveDetector("SlabMesh_SD", false));
    if (sdSlab) sdSlab->Save(tag, nev);
}
