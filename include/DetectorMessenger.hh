#pragma once
#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class DetectorMessenger : public G4UImessenger
{
public:
    DetectorMessenger(DetectorConstruction* det);
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand* cmd, G4String val) override;

private:
    DetectorConstruction*        fDet;
    G4UIdirectory*               fDetDir;
    G4UIcmdWithAString*          fMatCmd;
    G4UIcmdWithADoubleAndUnit*   fThickCmd;
};
