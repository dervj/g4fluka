# G4FlukaComparison — Geant4 analogue of g4fluka3.inp

## Overview

This Geant4 application replicates the FLUKA input file `g4fluka3.inp` for
benchmarking neutron transport through a slab target (50 cm Tungsten or
10.5 cm BPE).  It supports energy scans via macro files and writes CSV
output that directly mirrors FLUKA's USRBDX and USRYIELD scorer output.

---

## Build

```bash
mkdir build && cd build
cmake .. -DGeant4_DIR=/path/to/geant4/lib/Geant4-XX.Y.Z
make -j$(nproc)
```

Requires Geant4 ≥ 11.0 with the **G4NDL** neutron data library installed
(needed by `FTFP_BERT_HP` for thermal neutron transport).

---

## Run

**Single energy (interactive):**
```bash
./g4fluka
```

**Batch — Tungsten scan:**
```bash
./g4fluka scan_tungsten.mac
```

**Batch — BPE scan:**
```bash
./g4fluka scan_BPE.mac
```

Each run produces CSV files in the current directory (see Output section).

---

## FLUKA → Geant4 correspondence

| FLUKA card / concept          | Geant4 equivalent                                       |
|-------------------------------|---------------------------------------------------------|
| `DEFAULTS PRECISIO`           | `FTFP_BERT_HP` physics list (HP = high-precision data)  |
| `BEAM NEUTRON E=0.5 GeV`      | `G4ParticleGun` neutron, 500 MeV, `(0,0,+1)` direction |
| `BEAMPOS (0,0,-53)`           | `SetParticlePosition(0,0,-53*cm)`                       |
| `RPP SLAB` / `ASSIGNMA TUNGSTEN` | `G4Box("Slab")` + material `G4_W`                    |
| `MATERIAL BPE` / `COMPOUND`   | `G4Material("FlukaBPE")` with same mass fractions       |
| `MATERIAL POLYETHY`           | `G4Material("FlukaPolyethy")` 14.37% H + 85.63% C      |
| `MATERIAL BORON Z=5 ρ=2.37`   | `G4Material("FlukaBoron")` single-element               |
| `PART-THRES NEUTRON 1e-14 GeV`| Range cuts → ~0; HP library tracks neutrons to thermal  |
| `LOW-PWXS HYDROGEN h_polyethy`| G4NDL S(α,β) thermal scattering data (auto with HP)     |
| `USRBIN ENERGY` (EdepSlab)    | `MeshSD(EDEP)` 50×50×300 bins                           |
| `USRBIN NEUTRON` (NeutSlab)   | `MeshSD(NEUTRON_FLUX)` 100×100×105 bins                 |
| `USRBDX` BdxNeut/Prot/Phot…   | `SteppingAction` boundary crossing at SLAB→DSTRM        |
| `USRYIELD` YldNeut            | `SteppingAction::fYldNeut` d²N/dEdΩ, 18×100 bins       |
| `START 500000`                | `/run/beamOn 500000`                                    |

### Physics list notes

`FTFP_BERT_HP` is the closest Geant4 analogue to FLUKA's `PRECISIO` default
with neutron cross-section tables:

- **Bertini cascade** (0–10 GeV): comparable to FLUKA's PEANUT for hadronic
  inelastic reactions in the intermediate energy range.
- **HP data** (0–20 MeV): uses evaluated nuclear data files (ENDF/B-VIII,
  JEFF, JENDL) for elastic, inelastic, capture, and fission — analogous to
  FLUKA's low-energy neutron transport group structure (`LOW-PWXS`).
- S(α,β) thermal scattering: enabled automatically by G4NDL for H in
  polyethylene, corresponding to FLUKA's `LOW-PWXS HYDROGEN h_polyethy T=296`.
- EM physics: `G4EmStandardPhysics_option4` (same as FTFP_BERT_HP default),
  comparable to FLUKA `EMFCUT` thresholds.

### Known differences to be aware of

1. **Neutron group structure**: FLUKA uses a 260-group structure below 20 MeV;
   G4HP uses continuous cross sections from the same evaluated databases but
   with stochastic sampling.  Expect ≤5–10% differences in thermal flux.

2. **Hadronic models above 10 GeV**: not relevant here (beam ≤ 500 MeV) but
   note the transition from Bertini to FTFP at ~4 GeV.

3. **BPE S(α,β)**: FLUKA explicitly calls `h_polyethy` S(α,β) for H in
   polyethylene.  G4HP applies this automatically when the material contains H
   bonded in polyethylene — the `G4Material` name must contain "poly" or you
   must link the scattering table manually (see `G4NeutronHPThermalScattering`).

4. **Geometry**: FLUKA uses a `BLKBODY` sphere + `WORLD` RPP + VOID region.
   Geant4 uses a single `G4Box` world; the blackhole behaviour is implicit
   (particles leaving the world are killed).

5. **USRBDX normalisation**: FLUKA normalises to per primary per cm² per GeV
   (for Φ1 type).  The CSV output here uses the same normalisation (area =
   40 000 cm²).  Check that your FLUKA post-processing uses the same.

---

## Output files (per run tag)

| File                        | Content                                     | FLUKA analogue |
|-----------------------------|---------------------------------------------|----------------|
| `BdxNeut_<tag>.csv`         | Neutron spectrum at SLAB→DSTRM              | BdxNeut USRBDX |
| `BdxProt_<tag>.csv`         | Proton spectrum                             | BdxProt USRBDX |
| `BdxPhot_<tag>.csv`         | Photon spectrum                             | BdxPhot USRBDX |
| `BdxCHad_<tag>.csv`         | Charged hadron spectrum                     | BdxCHad USRBDX |
| `BdxElec_<tag>.csv`         | e+/e- spectrum                              | BdxElec USRBDX |
| `YldNeut_<tag>.csv`         | d²N/dEdΩ (18 θ × 100 E bins)               | YldNeut USRYIELD|
| `Edep_Zproj_<tag>.csv`      | Energy deposition Z-profile                 | EdepSlab USRBIN |
| `NeutFlux_Zproj_<tag>.csv`  | Neutron flux Z-profile                      | NeutSlab USRBIN |

---

## Plotting / comparison

```bash
# Compare G4 results with FLUKA for the 500 MeV tungsten case
python3 plot_comparison.py \
    --g4dir    ./build \
    --flukadir ./fluka_output \
    --tag      W50_500MeV \
    --fluka-tag W_500MeV \
    --outdir   plots/
```

FLUKA tab.lis files are produced in Flair by: selecting a USRBDX detector →
**Data** tab → **Export** → **tab-separated**.

---

## Extending the energy scan

Add additional `/beam/energy` + `/run/beamOn` lines in `scan_tungsten.mac`
or `scan_BPE.mac`.  The `G4FLUKA_TAG` environment variable controls the
output filename prefix.

To change THICK for BPE (e.g. 7 cm instead of 10.5 cm):
```
/det/setThickness 7 cm
```
