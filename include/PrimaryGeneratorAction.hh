#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class PrimaryMessenger;

// ─────────────────────────────────────────────────────────────────
// Pencil neutron beam along +Z, starting at z = -53 cm (FLUKA BEAMPOS)
// Energy set via /gun/energy or /beam/energy
// ─────────────────────────────────────────────────────────────────
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

    G4double GetBeamEnergy() const;
    void     SetBeamEnergy(G4double E);   // MeV

private:
    G4ParticleGun* fGun;
    PrimaryMessenger* fMessenger = nullptr;
};
