#!/usr/bin/env python3
"""
generate_all_plots.py
──────────────────────
Generate the full G4 plot set (USRBDX spectra, USRYIELD, USRBIN Z-profiles)
for every (material, energy) point produced by the HTCondor energy scan.

Scans condor/results/<TAG>/ for completed run directories (skipping any
"_mini" smoke-test directories) and calls the plotting functions already
defined in plot_comparison.py for each one, plus writes a summary_metrics.csv
of integrated quantities (leakage fluence per USRBDX particle, total
deposited energy, mean neutron flux) across all tags.

FLUKA output (in imports/) is named differently from the G4 tags -- e.g.
'BdxNeut_FLUKA_W_50cm_500MeV_proThOff_tab.lis' vs G4's
'BdxNeut_G4_W_50cm_500MeV.csv' -- so g4_tag_to_fluka_tag() below translates
between the two conventions. Per project decision, the '_proThOff' FLUKA
variant (PART-THRES card for protons deactivated) is used as canonical for
BOTH materials, since it's the only variant available for Tungsten.
USRYIELD is intentionally G4-only for now (no FLUKA overlay attempted) --
the available FLUKA yield exports don't match the G4 double-differential
format (see project notes).

Usage:
    python3 generate_all_plots.py
    python3 generate_all_plots.py --results-dir condor/results --outdir plots
    python3 generate_all_plots.py --flukadir imports
"""
import argparse
import csv
import glob
import os
import re

import plot_comparison as pc

NO_FLUKA_DIR = "/nonexistent"   # forces every FLUKA lookup to miss cleanly

_G4_TAG_RE = re.compile(r"(W|BPE)_([\d.]+)cm_([\d.]+)MeV$")

def g4_tag_to_fluka_tag(g4_tag, flukadir=None):
    """Map a G4 run tag (e.g. 'G4_W_50cm_500MeV', 'FlukaBPE_10.5cm_500MeV')
    to the corresponding FLUKA file tag (e.g. 'FLUKA_W_50cm_500MeV_proThOff',
    'FLUKA_BPE_10p5cm_500MeV_proThOff'). Returns None if the G4 tag doesn't
    match the expected material/thickness/energy pattern.

    '_proThOff' (PART-THRES for protons deactivated) is canonical, but at
    BPE/0.1MeV specifically the same setting was exported as
    '_protThresTest' instead (confirmed equivalent) -- if flukadir is given,
    probe for that fallback whenever the '_proThOff' file isn't present."""
    m = _G4_TAG_RE.search(g4_tag)
    if not m:
        return None
    material, thickness, energy = m.groups()
    thickness_p = thickness.replace(".", "p")
    energy_p = energy.replace(".", "p")
    base = f"FLUKA_{material}_{thickness_p}cm_{energy_p}MeV"
    candidates = [f"{base}_proThOff", f"{base}_protThresTest"]

    if flukadir:
        for cand in candidates:
            if os.path.exists(os.path.join(flukadir, f"BdxNeut_{cand}_tab.lis")):
                return cand

    return candidates[0]


def find_tag(result_dir):
    """Derive the run tag (e.g. G4_W_50cm_500MeV) from an Edep_Zproj_*.csv
    filename already present in the result directory."""
    matches = glob.glob(os.path.join(result_dir, "Edep_Zproj_*.csv"))
    if not matches:
        return None
    fname = os.path.basename(matches[0])
    return fname[len("Edep_Zproj_"):-len(".csv")]


_BEAMON_RE = re.compile(r"/run/beamOn\s+(\d+)")

def get_n_primaries(run_dir_name, macros_dir="macros/energy_points"):
    """Look up the number of primaries for a run from its macro file --
    the result directory name (e.g. 'W50_500MeV') matches the macro
    basename by construction. Needed to reconstruct per-bin Poisson counts
    for USRBDX statistical error bars (the CSVs only store the normalized
    value, but the raw count is exactly recoverable given N)."""
    macro_path = os.path.join(macros_dir, f"{run_dir_name}.mac")
    if not os.path.exists(macro_path):
        return None
    with open(macro_path) as f:
        m = _BEAMON_RE.search(f.read())
    return int(m.group(1)) if m else None


def generate_for_tag(g4dir, tag, outdir, flukadir):
    print(f"--- {tag} ---")
    fluka_tag = g4_tag_to_fluka_tag(tag, flukadir) if flukadir != NO_FLUKA_DIR else tag
    if flukadir != NO_FLUKA_DIR and fluka_tag is None:
        print(f"  [warn] could not derive a FLUKA tag from '{tag}'; skipping FLUKA overlay")
        fluka_tag = tag

    n_primaries = get_n_primaries(os.path.basename(g4dir))
    if n_primaries is None:
        print(f"  [warn] could not determine N primaries for '{tag}'; no G4 error bars")

    for pname in pc.PARTICLE_LABELS:
        pc.compare_one(pname, g4dir, tag, flukadir, fluka_tag, outdir, n_primaries)

    pc.plot_yield(g4dir, tag, outdir)

    for qname in pc.USRBIN_LABELS:
        pc.plot_usrbin(qname, g4dir, tag, flukadir, fluka_tag, outdir)

    return pc.compute_summary_metrics(tag, g4dir, flukadir, fluka_tag)


def write_summary_csv(rows, path):
    fieldnames = []
    for row in rows:
        for k in row:
            if k not in fieldnames:
                fieldnames.append(k)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="condor/results",
                     help="Directory containing one subdirectory per run tag")
    ap.add_argument("--outdir", default="plots",
                     help="Output directory for PNGs (one subfolder per tag)")
    ap.add_argument("--flukadir", default=None,
                     help="Directory with FLUKA ASCII output, named with the "
                          "same tag as the G4 CSVs. Omit to skip FLUKA overlay.")
    args = ap.parse_args()

    flukadir = args.flukadir if args.flukadir else NO_FLUKA_DIR

    run_dirs = sorted(
        d for d in glob.glob(os.path.join(args.results_dir, "*"))
        if os.path.isdir(d) and not os.path.basename(d).endswith("_mini")
    )
    if not run_dirs:
        print(f"No run directories found under {args.results_dir}")
        return

    summary_rows = []
    for d in run_dirs:
        tag = find_tag(d)
        if tag is None:
            print(f"  [skip] {d}: no Edep_Zproj_*.csv found")
            continue
        outdir = os.path.join(args.outdir, os.path.basename(d))
        os.makedirs(outdir, exist_ok=True)
        summary_rows.append(generate_for_tag(d, tag, outdir, flukadir))

    write_summary_csv(summary_rows, os.path.join(args.outdir, "summary_metrics.csv"))
    print(f"\nDone. Plots written under {args.outdir}/<run_tag>/")


if __name__ == "__main__":
    main()
