#pragma once
#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

class PrimaryMessenger : public G4UImessenger
{
public:
    PrimaryMessenger(PrimaryGeneratorAction* gun);
    ~PrimaryMessenger() override;
    void SetNewValue(G4UIcommand* cmd, G4String val) override;

private:
    PrimaryGeneratorAction*    fGun;
    G4UIdirectory*             fBeamDir;
    G4UIcmdWithADoubleAndUnit* fEnergyCmd;
};
