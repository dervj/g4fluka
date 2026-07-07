#!/usr/bin/env python3
"""
plot_comparison.py
──────────────────
Plot G4 boundary-crossing spectra and compare with FLUKA USRBDX output.

Usage:
    python3 plot_comparison.py --g4dir ./results --flukadir ./fluka_results \
            --tag W50_500MeV --fluka-tag W_500MeV

The script reads:
  G4:    BdxNeut_<tag>.csv, BdxProt_<tag>.csv, etc.
  FLUKA: BdxNeut_<fluka-tag>_tab.lis  (from Flair → Data → Export)
         or any CSV with columns Elow,Ehigh,value

Outputs: comparison plots as PNG.
"""
import argparse
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 11,
    "axes.grid": True,
    "grid.alpha": 0.3,
})

PARTICLE_LABELS = {
    "BdxNeut": ("Neutron",       "tab:blue"),
    "BdxProt": ("Proton",        "tab:orange"),
    "BdxPhot": ("Photon",        "tab:green"),
    "BdxCHad": ("Charged hadron","tab:red"),
    "BdxElec": ("e+/e-",         "tab:purple"),
}

def load_csv(path):
    """Load a CSV with columns: Elow, Ehigh, value (skips # lines)."""
    data = np.loadtxt(path, delimiter=",", comments="#")
    return data[:,0], data[:,1], data[:,2]

_N_INTERVALS_RE = re.compile(r"N\.\s*of\s+(?:energy|x1)\s+intervals\s+(\d+)", re.IGNORECASE)

def load_fluka_lis(path):
    """
    Load a FLUKA USRBDX/USRYIELD ASCII export (Flair "tab.lis" or the
    equivalent rfluka post-processing dump).

    These files are NOT simply "every numeric line is a data row": after the
    primary "integrated over solid angle" block (N rows, one per energy/x1
    bin) they often continue with a SECOND block for each solid-angle
    interval, repeating the same bin edges with a different column count.
    Blindly parsing every >=3-column numeric line would silently concatenate
    both blocks and corrupt the bin sequence (duplicate/out-of-order x
    values). So: find the declared interval count from the header comment
    ("# N. of energy intervals NNN" / "# N. of x1 intervals NNN") and read
    exactly that many data rows following it, ignoring anything after.
    """
    with open(path) as f:
        lines = f.readlines()

    n_expected = None
    start_idx = 0
    for i, line in enumerate(lines):
        m = _N_INTERVALS_RE.search(line)
        if m:
            n_expected = int(m.group(1))
            start_idx = i + 1
            break

    rows = []
    for line in lines[start_idx:]:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("*"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                continue
        if n_expected is not None and len(rows) >= n_expected:
            break

    data = np.array(rows)
    return data[:,0], data[:,1], data[:,2]

def plot_spectrum(ax, elow, ehigh, vals, label, color, ls="-"):
    emid = np.sqrt(elow * ehigh)   # geometric midpoint for log-scale
    ax.stairs(vals, np.append(elow, ehigh[-1]),
              label=label, color=color, linestyle=ls, linewidth=1.5)

def rebin_loglog(x_target, x_src, y_src):
    """Log-log interpolate y_src(x_src) onto x_target. Returns NaN outside
    the source's range (no extrapolation) or where source values are <=0."""
    mask = y_src > 0
    if mask.sum() < 2:
        return np.full_like(x_target, np.nan)
    logx_src = np.log10(x_src[mask])
    logy_src = np.log10(y_src[mask])
    order = np.argsort(logx_src)
    logx_src, logy_src = logx_src[order], logy_src[order]
    logx_t = np.log10(x_target)
    logy_t = np.interp(logx_t, logx_src, logy_src, left=np.nan, right=np.nan)
    return 10**logy_t

def rebin_linear(x_target, x_src, y_src):
    """Linear interpolate y_src(x_src) onto x_target; NaN outside range."""
    return np.interp(x_target, x_src, y_src, left=np.nan, right=np.nan)

def make_ratio_figure():
    """Two-panel figure: overlay on top, G4/FLUKA ratio beneath."""
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7, 6.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})
    return fig, ax, axr

def compare_one(pname, g4dir, tag, flukadir, fluka_tag, outdir):
    label, color = PARTICLE_LABELS.get(pname, (pname, "black"))

    g4_file = os.path.join(g4dir, f"{pname}_{tag}.csv")
    if not os.path.exists(g4_file):
        print(f"  [skip] {g4_file} not found")
        return

    el_g4, eh_g4, val_g4 = load_csv(g4_file)
    emid_g4 = np.sqrt(el_g4 * eh_g4)

    # Try to load FLUKA result
    fl_data = None
    for ext in ["_tab.lis", ".csv", ".lis"]:
        fl_file = os.path.join(flukadir, f"{pname}_{fluka_tag}{ext}")
        if os.path.exists(fl_file):
            try:
                if ext == ".csv":
                    fl_data = load_csv(fl_file)
                else:
                    fl_data = load_fluka_lis(fl_file)
            except Exception as e:
                print(f"  [warn] Could not load FLUKA file {fl_file}: {e}")
            break

    if fl_data is not None:
        el_fl, eh_fl, val_fl = fl_data
        emid_fl = np.sqrt(el_fl * eh_fl)
        val_fl_on_g4 = rebin_loglog(emid_g4, emid_fl, val_fl)
        ratio = val_g4 / val_fl_on_g4

        fig, ax, axr = make_ratio_figure()
        plot_spectrum(ax, el_g4, eh_g4, val_g4, f"Geant4 ({tag})", color, "-")
        plot_spectrum(ax, el_fl, eh_fl, val_fl, f"FLUKA ({fluka_tag})", "black", "--")
        axr.plot(emid_g4, ratio, color=color, marker=".", markersize=3, linewidth=0.8)
        axr.axhline(1.0, color="black", linewidth=0.8, linestyle=":")
        axr.set_ylim(0, 2)
        axr.set_ylabel("G4 / FLUKA")
        axr.set_xlabel("Kinetic Energy [GeV]")
        ax.set_ylabel(r"$d\Phi/dE$ [cm$^{-2}$ GeV$^{-1}$ primary$^{-1}$]")
        ax.set_title(f"{label} spectrum – SLAB→DSTRM boundary\nTag: {tag}")
        ax.legend()
        ax.set_ylim(bottom=1e-12)
        ax.set_xscale("log")
        ax.set_yscale("log")
    else:
        fig, ax = plt.subplots(figsize=(7,5))
        plot_spectrum(ax, el_g4, eh_g4, val_g4, f"Geant4 ({tag})", color, "-")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Kinetic Energy [GeV]")
        ax.set_ylabel(r"$d\Phi/dE$ [cm$^{-2}$ GeV$^{-1}$ primary$^{-1}$]")
        ax.set_title(f"{label} spectrum – SLAB→DSTRM boundary\nTag: {tag}")
        ax.legend()
        ax.set_ylim(bottom=1e-12)

    out = os.path.join(outdir, f"{pname}_{tag}_comparison.png")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")

def plot_yield(g4dir, tag, outdir):
    fname = os.path.join(g4dir, f"YldNeut_{tag}.csv")
    if not os.path.exists(fname):
        print(f"  [skip] {fname} not found")
        return

    data = np.loadtxt(fname, delimiter=",", comments="#")
    # columns: theta_lo, theta_hi, Elo, Ehi, value
    theta_vals = np.unique(data[:,0])
    e_lo       = np.unique(data[:,2])

    fig, ax = plt.subplots(figsize=(8,6))
    cmap = plt.cm.viridis
    for i, th in enumerate(theta_vals):
        mask = data[:,0] == th
        chunk = data[mask]
        emid = np.sqrt(chunk[:,2]*chunk[:,3])
        ax.plot(emid, chunk[:,4],
                color=cmap(i/len(theta_vals)),
                label=f"θ={th:.0f}–{data[mask,1][0]:.0f}°",
                linewidth=1.)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Kinetic Energy [GeV]")
    ax.set_ylabel(r"$d^2N/dEd\Omega$ [GeV$^{-1}$ sr$^{-1}$ primary$^{-1}$]")
    ax.set_title(f"Neutron double-differential yield – {tag}")
    ax.legend(fontsize=7, ncol=3, loc="upper right")
    ax.set_ylim(bottom=1e-14)

    out = os.path.join(outdir, f"YldNeut_{tag}.png")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")

USRBIN_LABELS = {
    "Edep_Zproj":     ("Energy deposition", "tab:red",  r"Energy deposit [GeV primary$^{-1}$]"),
    "NeutFlux_Zproj": ("Neutron flux",      "tab:blue", r"Track-length density [cm primary$^{-1}$ cm$^{-1}_{\rm depth}$]"),
}

# FLUKA's USRBIN export uses different detector names than our G4 CSVs, and
# a different per-bin convention (a Z-linear density, not a per-bin total):
#   EdepSlab_*_Zprof.lis:  GeV/cm/primary  -> multiply by dz to get a per-bin
#                          total, matching G4's Edep_Zproj convention.
#   NeutSlab_*_Zprof.lis:  cm/primary/cm_depth, i.e. FLUKA already sums
#                          (density * voxel_area) over the XY plane. Our G4
#                          MeshSD Z-projection instead sums raw per-voxel
#                          densities (no area weighting), so it must be
#                          multiplied by the (fixed, material-independent)
#                          NeutFlux mesh voxel area to match FLUKA's
#                          convention. See DetectorConstruction.cc: the
#                          NeutFlux/NeutSlab grid is always 100x100 bins
#                          over [-100,100] cm in both x and y => dx=dy=2cm.
FLUKA_USRBIN_PREFIX = {"Edep_Zproj": "EdepSlab", "NeutFlux_Zproj": "NeutSlab"}
NEUT_VOXEL_AREA_CM2 = 2.0 * 2.0

def plot_usrbin(qname, g4dir, tag, flukadir, fluka_tag, outdir):
    """Plot the Z-profile (XY-integrated) of a USRBIN-equivalent mesh,
    e.g. Edep_Zproj_<tag>.csv or NeutFlux_Zproj_<tag>.csv."""
    label, color, ylabel = USRBIN_LABELS[qname]

    g4_file = os.path.join(g4dir, f"{qname}_{tag}.csv")
    if not os.path.exists(g4_file):
        print(f"  [skip] {g4_file} not found")
        return

    zlo_g4, zhi_g4, val_g4 = load_csv(g4_file)
    if qname == "NeutFlux_Zproj":
        val_g4 = val_g4 * NEUT_VOXEL_AREA_CM2
    zmid_g4 = 0.5*(zlo_g4 + zhi_g4)

    # Try to load FLUKA USRBIN result (different prefix/extension than Bdx)
    fluka_prefix = FLUKA_USRBIN_PREFIX.get(qname, qname)
    fl_data = None
    for ext in ["_Zprof.lis", "_tab.lis", ".csv", ".lis"]:
        fl_file = os.path.join(flukadir, f"{fluka_prefix}_{fluka_tag}{ext}")
        if os.path.exists(fl_file):
            try:
                if ext == ".csv":
                    fl_data = load_csv(fl_file)
                else:
                    fl_data = load_fluka_lis(fl_file)
            except Exception as e:
                print(f"  [warn] Could not load FLUKA file {fl_file}: {e}")
            break

    if fl_data is not None:
        zlo_fl, zhi_fl, val_fl = fl_data
        if qname == "Edep_Zproj":
            val_fl = val_fl * (zhi_fl - zlo_fl)
        zmid_fl = 0.5*(zlo_fl + zhi_fl)
        val_fl_on_g4 = rebin_linear(zmid_g4, zmid_fl, val_fl)
        ratio = val_g4 / val_fl_on_g4

        fig, ax, axr = make_ratio_figure()
        ax.stairs(val_g4, np.append(zlo_g4, zhi_g4[-1]),
                  label=f"Geant4 ({tag})", color=color, linewidth=1.5)
        ax.stairs(val_fl, np.append(zlo_fl, zhi_fl[-1]),
                  label=f"FLUKA ({fluka_tag})", color="black",
                  linestyle="--", linewidth=1.5)
        axr.plot(zmid_g4, ratio, color=color, marker=".", markersize=3, linewidth=0.8)
        axr.axhline(1.0, color="black", linewidth=0.8, linestyle=":")
        axr.set_ylim(0, 2)
        axr.set_ylabel("G4 / FLUKA")
        axr.set_xlabel("z [cm]  (slab entrance = 0)")
        ax.set_yscale("log")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{label} — Z profile\nTag: {tag}")
        ax.legend()
    else:
        fig, ax = plt.subplots(figsize=(7,5))
        ax.stairs(val_g4, np.append(zlo_g4, zhi_g4[-1]),
                  label=f"Geant4 ({tag})", color=color, linewidth=1.5)
        ax.set_yscale("log")
        ax.set_xlabel("z [cm]  (slab entrance = 0)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{label} — Z profile\nTag: {tag}")
        ax.legend()

    out = os.path.join(outdir, f"{qname}_{tag}_comparison.png")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")

def _load_any(prefix, tag, directory):
    """Try each supported extension for prefix_tag.*; return (elo,ehi,val)
    or None if nothing found/loadable."""
    for ext in ["_Zprof.lis", "_tab.lis", ".csv", ".lis"]:
        path = os.path.join(directory, f"{prefix}_{tag}{ext}")
        if os.path.exists(path):
            try:
                return load_csv(path) if ext == ".csv" else load_fluka_lis(path)
            except Exception as e:
                print(f"  [warn] Could not load {path}: {e}")
                return None
    return None

def compute_summary_metrics(tag, g4dir, flukadir=None, fluka_tag=None):
    """Integrated/summary metrics for one run tag: total leakage fluence per
    USRBDX particle [cm^-2 primary^-1], total deposited energy [GeV primary^-1],
    and mean neutron track-length flux [cm^-2 primary^-1]. Includes FLUKA
    values and the G4/FLUKA ratio wherever a matching FLUKA file is found."""
    row = {"tag": tag}

    for pname in PARTICLE_LABELS:
        g4_data = _load_any(pname, tag, g4dir)
        if g4_data is None:
            continue
        elo, ehi, val = g4_data
        g4_integral = float(np.sum(val * (ehi - elo)))
        row[f"{pname}_G4_integral_fluence"] = g4_integral

        if flukadir and fluka_tag:
            fl_data = _load_any(pname, fluka_tag, flukadir)
            if fl_data is not None:
                fl_elo, fl_ehi, fl_val = fl_data
                fl_integral = float(np.sum(fl_val * (fl_ehi - fl_elo)))
                row[f"{pname}_FLUKA_integral_fluence"] = fl_integral
                row[f"{pname}_ratio_G4_over_FLUKA"] = (
                    g4_integral / fl_integral if fl_integral else float("nan"))

    edep_data = _load_any("Edep_Zproj", tag, g4dir)
    if edep_data is not None:
        _, _, val = edep_data
        row["Edep_total_GeV_per_primary"] = float(np.sum(val))
        if flukadir and fluka_tag:
            fl_data = _load_any(FLUKA_USRBIN_PREFIX["Edep_Zproj"], fluka_tag, flukadir)
            if fl_data is not None:
                fl_zlo, fl_zhi, fl_val = fl_data
                fl_total = float(np.sum(fl_val * (fl_zhi - fl_zlo)))  # density -> per-bin total
                row["Edep_FLUKA_total_GeV_per_primary"] = fl_total
                row["Edep_ratio_G4_over_FLUKA"] = (
                    row["Edep_total_GeV_per_primary"] / fl_total if fl_total else float("nan"))

    neut_data = _load_any("NeutFlux_Zproj", tag, g4dir)
    if neut_data is not None:
        _, _, val = neut_data
        val = val * NEUT_VOXEL_AREA_CM2  # match FLUKA's area-integrated convention
        row["NeutFlux_mean_per_primary"] = float(np.mean(val))
        if flukadir and fluka_tag:
            fl_data = _load_any(FLUKA_USRBIN_PREFIX["NeutFlux_Zproj"], fluka_tag, flukadir)
            if fl_data is not None:
                _, _, fl_val = fl_data
                fl_mean = float(np.mean(fl_val))
                row["NeutFlux_FLUKA_mean_per_primary"] = fl_mean
                row["NeutFlux_ratio_G4_over_FLUKA"] = (
                    row["NeutFlux_mean_per_primary"] / fl_mean if fl_mean else float("nan"))

    return row

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g4dir",     default=".", help="Directory with G4 CSV files")
    ap.add_argument("--flukadir",  default=".", help="Directory with FLUKA output files")
    ap.add_argument("--tag",       default="W50_500MeV", help="G4 run tag")
    ap.add_argument("--fluka-tag", default="W_500MeV",   help="FLUKA run tag")
    ap.add_argument("--outdir",    default="plots",      help="Output directory for PNGs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"G4 dir   : {args.g4dir}")
    print(f"FLUKA dir: {args.flukadir}")
    print(f"Tag      : {args.tag}  vs  {args.fluka_tag}")

    for pname in PARTICLE_LABELS:
        compare_one(pname, args.g4dir, args.tag,
                    args.flukadir, args.fluka_tag, args.outdir)

    plot_yield(args.g4dir, args.tag, args.outdir)

    for qname in USRBIN_LABELS:
        plot_usrbin(qname, args.g4dir, args.tag,
                    args.flukadir, args.fluka_tag, args.outdir)

if __name__ == "__main__":
    main()
