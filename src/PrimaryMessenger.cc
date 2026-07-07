#include "PrimaryMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

PrimaryMessenger::PrimaryMessenger(PrimaryGeneratorAction* gun) : fGun(gun)
{
    fBeamDir = new G4UIdirectory("/beam/");
    fBeamDir->SetGuidance("Beam configuration");

    fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/beam/energy", this);
    fEnergyCmd->SetGuidance("Set neutron beam kinetic energy");
    fEnergyCmd->SetParameterName("E", false);
    fEnergyCmd->SetDefaultUnit("MeV");
    fEnergyCmd->SetUnitCandidates("eV keV MeV GeV");
    fEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryMessenger::~PrimaryMessenger()
{
    delete fEnergyCmd;
    delete fBeamDir;
}

void PrimaryMessenger::SetNewValue(G4UIcommand* cmd, G4String val)
{
    if (cmd == fEnergyCmd)
        fGun->SetBeamEnergy(fEnergyCmd->GetNewDoubleValue(val));
}
