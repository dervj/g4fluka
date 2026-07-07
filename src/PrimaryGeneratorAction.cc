#include "PrimaryGeneratorAction.hh"
#include "PrimaryMessenger.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// ─────────────────────────────────────────────────────────────────
// FLUKA BEAM card: 0.5 GeV neutron, pencil beam along +Z
// FLUKA BEAMPOS:   position (0,0,-53) cm
// ─────────────────────────────────────────────────────────────────

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fGun = new G4ParticleGun(1);   // 1 particle per event

    G4ParticleDefinition* neutron =
        G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    fGun->SetParticleDefinition(neutron);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fGun->SetParticlePosition(G4ThreeVector(0., 0., -53.*cm));
    fGun->SetParticleEnergy(500.*MeV);   // FLUKA: E=0.5 GeV kinetic

    fMessenger = new PrimaryMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
    delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    fGun->GeneratePrimaryVertex(event);
}

G4double PrimaryGeneratorAction::GetBeamEnergy() const
{
    return fGun->GetParticleEnergy();
}

void PrimaryGeneratorAction::SetBeamEnergy(G4double E)
{
    fGun->SetParticleEnergy(E);
}
