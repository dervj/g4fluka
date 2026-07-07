#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
    : fDet(det)
{
    fDetDir = new G4UIdirectory("/det/");
    fDetDir->SetGuidance("Detector configuration commands");

    fMatCmd = new G4UIcmdWithAString("/det/setMaterial", this);
    fMatCmd->SetGuidance("Set slab material: Tungsten or BPE");
    fMatCmd->SetParameterName("material", false);
    fMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fThickCmd = new G4UIcmdWithADoubleAndUnit("/det/setThickness", this);
    fThickCmd->SetGuidance("Set slab thickness (e.g. 50 cm or 10.5 cm)");
    fThickCmd->SetParameterName("thickness", false);
    fThickCmd->SetDefaultUnit("cm");
    fThickCmd->SetUnitCandidates("mm cm m");
    fThickCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
    delete fMatCmd;
    delete fThickCmd;
    delete fDetDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* cmd, G4String val)
{
    if (cmd == fMatCmd)
        fDet->SetSlabMaterial(val);
    else if (cmd == fThickCmd)
        fDet->SetSlabThickness(fThickCmd->GetNewDoubleValue(val));
}
