// ─────────────────────────────────────────────────────────────────
// g4fluka.cc  –  Main driver
//
// Geant4 analogue of the FLUKA g4fluka3.inp setup.
// Usage:
//   ./g4fluka run_W_500MeV.mac       (batch)
//   ./g4fluka                         (interactive)
// ─────────────────────────────────────────────────────────────────
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

// Physics list – FTFP_BERT_HP is the standard for neutron work
// (HP = high-precision thermal neutrons; analogous to FLUKA PRECISIO + LOW-PWXS)
#include "FTFP_BERT_HP.hh"
#include "G4HadronicParameters.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"

int main(int argc, char** argv)
{
    // ── Run manager ───────────────────────────────────────────────
    auto* runMgr = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::SerialOnly);

    // ── Detector ─────────────────────────────────────────────────
    auto* det = new DetectorConstruction();
    runMgr->SetUserInitialization(det);

    // ── Physics list ──────────────────────────────────────────────
    // FTFP_BERT_HP:
    //  • BERT  = Bertini cascade (< ~10 GeV)  ≈ FLUKA intra-nuclear cascade
    //  • FTFP  = Fritiof string model (> ~4 GeV)
    //  • HP    = G4NDL high-precision neutron data (thermal–20 MeV)
    //            analogous to FLUKA LOW-PWXS + h_polyethy S(α,β)
    // Thresholds: production cut set globally below to match FLUKA
    //   PART-THRES NEUTRON 1e-14 GeV → transport down to 0 kinetic energy
    //   EMFCUT e+/e- 1 MeV, gamma 0.1 MeV  (FLUKA EMFCUT cards)
    auto* physicsList = new FTFP_BERT_HP();
    physicsList->SetDefaultCutValue(0.1*mm);
    runMgr->SetUserInitialization(physicsList);

    // ── User actions ─────────────────────────────────────────────
    auto* primaryGen = new PrimaryGeneratorAction();
    runMgr->SetUserAction(primaryGen);

    auto* runAction = new RunAction();
    runMgr->SetUserAction(runAction);

    auto* eventAction = new EventAction();
    runMgr->SetUserAction(eventAction);

    // SteppingAction needs the detector for geometry queries
    auto* steppingAction = new SteppingAction(det);
    runMgr->SetUserAction(steppingAction);

    // ── UI/Vis or batch ──────────────────────────────────────────
    G4UImanager* UI = G4UImanager::GetUIpointer();

    if (argc > 1) {
        // Batch mode: execute macro
        G4String macFile = argv[1];
        UI->ApplyCommand("/control/execute " + macFile);
    } else {
        // Interactive mode
        G4UIExecutive* uiExec = new G4UIExecutive(argc, argv);
        G4VisManager* visMgr = new G4VisExecutive();
        visMgr->Initialize();
        UI->ApplyCommand("/control/execute macros/init_vis.mac");
        uiExec->SessionStart();
        delete uiExec;
        delete visMgr;
    }

    delete runMgr;
    return 0;
}
