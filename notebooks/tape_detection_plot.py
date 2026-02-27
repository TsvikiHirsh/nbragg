"""
Tape Detection Demo
===================
Shows that ~100 µm painter's tape (cellulose) on a pepper can is detectable
in neutron TOF transmission by comparing fits with and without hydrous material.

Strategy
--------
• Fit the pure iron reference to get the "correct" iron cross-section parameters.
• Create a synthetic "pepper can + tape" dataset by multiplying the noiseless
  iron model by the cellulose tape attenuation and adding realistic noise.
• Panel 1: apply iron-only parameters (from reference fit) to the tape dataset
            → systematic residuals because cellulose is unaccounted for.
• Panel 2: fit only the cellulose thickness (all iron params fixed)
            → residuals become flat; fitted thickness ≈ 100 µm.
"""
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import lmfit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import NCrystal as NC
import nbragg
from scipy.ndimage import uniform_filter1d   # for smoothed residual trend

RNG = np.random.default_rng(42)

# ─────────────────────────────────────────────
# 1.  Materials
# ─────────────────────────────────────────────
nbragg.register_material("Fe_sg229_Iron-alpha_CrysExtn1.ncmat")
nbragg.register_material("Cellulose_C6O5H10.ncmat")

# ─────────────────────────────────────────────
# 2.  Iron reference fit (no background, no extinction)
# ─────────────────────────────────────────────
_raw = pd.read_csv("iron_powder_transmission.csv", sep=r"\s+",
                   header=0, names=["wl", "trans", "err"])
wl  = _raw["wl"].values
err = _raw["err"].values

xs_iron  = nbragg.CrossSection(iron=nbragg.materials["Fe_sg229_Iron-alpha_CrysExtn1.ncmat"])
model_fe = nbragg.TransmissionModel(xs_iron, vary_background=False, vary_extinction=False)

data_pure = nbragg.Data.from_transmission(_raw.rename(columns={"wl": "wavelength"}))
print("Fitting pure iron reference …")
res_pure   = model_fe.fit(data_pure, wlmin=0.5, wlmax=9.0)
print(f"  Iron reference reduced χ² = {res_pure.redchi:.1f}")

params_ref   = res_pure.params.valuesdict()      # fixed "truth" iron params
T_iron_ref   = model_fe.transmission(wl, **params_ref)   # noiseless iron curve

# ─────────────────────────────────────────────
# 3.  Cellulose (painter's tape)
# ─────────────────────────────────────────────
cel_mat      = NC.load("Cellulose_C6O5H10.ncmat")
n_cel        = cel_mat.info.factor_macroscopic_xs
cel_xs       = (cel_mat.scatter.xsect(wl=wl, direction=(0,0,1)) +
                cel_mat.absorption.xsect(wl=wl, direction=(0,0,1)))

TAPE_TRUE_CM = 0.010            # 100 µm
T_tape_true  = np.exp(-cel_xs * n_cel * TAPE_TRUE_CM)
print(f"\nTape effect: {(1-T_tape_true.min())*100:.1f}% at λ={wl[np.argmin(T_tape_true)]:.1f} Å")

# ─────────────────────────────────────────────
# 4.  Synthetic tape dataset
# ─────────────────────────────────────────────
T_data   = T_iron_ref * T_tape_true + RNG.normal(0., err)
err_data = err.copy()

mask     = (wl >= 0.5) & (wl <= 9.0)
wl_fit   = wl[mask]
T_fit    = T_data[mask]
err_fit  = err_data[mask]
cel_fit  = cel_xs[mask]

# ─────────────────────────────────────────────
# 5.  Panel 1 model — iron reference (no fitting)
#     Iron params are FIXED at reference values.
#     Residuals reveal the cellulose signature.
# ─────────────────────────────────────────────
T_panel1   = T_iron_ref
resid1_raw = (T_data - T_panel1) / err_data
chi2_iron  = float(np.mean(resid1_raw[mask]**2))
print(f"\nIron-reference (fixed) reduced χ² on tape data = {chi2_iron:.1f}")

# ─────────────────────────────────────────────
# 6.  Panel 2 — fit only cellulose thickness
#     Iron params stay fixed; only t_cel varies.
# ─────────────────────────────────────────────
params_cel = lmfit.Parameters()
params_cel.add("t_cel", value=0.005, min=1e-5, max=0.5, vary=True)  # start at 50 µm

def residual_cel_only(p):
    T_ce = np.exp(-cel_fit * n_cel * p["t_cel"].value)
    return (T_fit - T_panel1[mask] * T_ce) / err_fit

mini_cel  = lmfit.Minimizer(residual_cel_only, params_cel)
res_cel   = mini_cel.minimize(method="leastsq")
t_cel_fit = res_cel.params["t_cel"].value
chi2_comb = res_cel.redchi
print(f"Iron + cellulose (fit t_cel only) reduced χ² = {chi2_comb:.2f}")
print(f"Fitted tape thickness = {t_cel_fit*1e4:.0f} µm  (true = {TAPE_TRUE_CM*1e4:.0f} µm)")

T_panel2   = T_panel1 * np.exp(-cel_xs * n_cel * t_cel_fit)
resid2_raw = (T_data - T_panel2) / err_data

# Smoothed residual trend (200-pt rolling mean to reveal systematic bias)
# Use 'nearest' edge mode to avoid edge artefacts
SMOOTH = 200
trend1 = uniform_filter1d(resid1_raw, SMOOTH, mode="nearest")
trend2 = uniform_filter1d(resid2_raw, SMOOTH, mode="nearest")
# Also compute the THEORETICAL cellulose residual (what the fit should recover)
theory_resid = -(1 - T_tape_true) * T_iron_ref / err_data   # ≈ negative, deepens with λ

# ─────────────────────────────────────────────
# 7.  Figure
# ─────────────────────────────────────────────
RED    = "#c0392b"
GREEN  = "#27ae60"
GREY   = "0.38"
WLMIN, WLMAX = 0.5, 9.0

fig = plt.figure(figsize=(14, 9))
fig.patch.set_facecolor("white")
gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1], hspace=0.07, wspace=0.30,
                       left=0.07, right=0.97, top=0.88, bottom=0.09)
ax_l  = fig.add_subplot(gs[0, 0])
ax_lr = fig.add_subplot(gs[1, 0], sharex=ax_l)
ax_r  = fig.add_subplot(gs[0, 1])
ax_rr = fig.add_subplot(gs[1, 1], sharex=ax_r)

step = max(1, len(wl) // 500)

def draw_panel(ax_main, ax_res, T_model, resid_raw, trend,
               label, color, title, chi2_label, annotation=None):

    # ── Main plot ──
    ax_main.errorbar(wl[::step], T_data[::step], yerr=err_data[::step],
                     fmt=".", color=GREY, ms=2.5, lw=0.5, alpha=0.5,
                     label="Data  (iron can + tape)", zorder=1)
    ax_main.plot(wl, T_model, color=color, lw=2.2, label=label, zorder=3)

    ax_main.set_xlim(WLMIN, WLMAX)
    ax_main.set_ylabel("Neutron Transmission", fontsize=10)
    ax_main.legend(fontsize=8.5, loc="upper right")
    ax_main.tick_params(labelbottom=False)
    ax_main.set_title(title, fontsize=11, fontweight="bold", pad=5)
    ax_main.text(0.03, 0.93, chi2_label,
                 transform=ax_main.transAxes, fontsize=9.5, va="top",
                 bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="0.6"))
    if annotation:
        ax_main.text(0.97, 0.08, annotation,
                     transform=ax_main.transAxes, ha="right", va="bottom",
                     fontsize=9.5, color=color,
                     bbox=dict(boxstyle="round,pad=0.3", fc="white",
                               ec=color, alpha=0.9))

    # ── Residuals ──
    in_range = (wl >= WLMIN) & (wl <= WLMAX)
    ax_res.axhline(0, color="0.3", lw=0.8, ls="--")
    ax_res.fill_between(wl, -1, 1, color="lightgreen", alpha=0.35)
    ax_res.plot(wl, resid_raw, ".", color=color, ms=1.0, alpha=0.25)
    ax_res.plot(wl, trend,     color=color,  lw=2.2,  alpha=0.95,
                label="Smoothed trend")
    ax_res.legend(fontsize=7.5, loc="lower left")


    ax_res.set_ylabel("Residual (σ)", fontsize=9)
    ax_res.set_xlabel("Wavelength (Å)", fontsize=10)
    ax_res.set_xlim(WLMIN, WLMAX)
    ylim = max(4, np.percentile(np.abs(resid_raw[in_range]), 99) * 0.7)
    ax_res.set_ylim(-ylim, ylim)
    ax_res.yaxis.set_major_locator(plt.MultipleLocator(2))


draw_panel(ax_l, ax_lr,
           T_model=T_panel1,
           resid_raw=resid1_raw,
           trend=trend1,
           label="Iron-only prediction",
           color=RED,
           title="Without hydrous material — poor fit",
           chi2_label=f"χ² / d.o.f. = {chi2_iron:.1f}  (expect 1.0 for good fit)")

# Overlay theoretical cellulose residual on the left residual panel
ax_lr.plot(wl, theory_resid, color="darkorange", lw=1.5, ls="--", alpha=0.85,
           label=f"Expected H-scattering signal\n({TAPE_TRUE_CM*1e4:.0f} µm cellulose)")
ax_lr.legend(fontsize=7, loc="lower left", ncol=2)

tape_note = (f"Fitted tape:  {t_cel_fit*1e4:.0f} µm\n"
             f"True value:  {TAPE_TRUE_CM*1e4:.0f} µm  (painter's tape)")
draw_panel(ax_r, ax_rr,
           T_model=T_panel2,
           resid_raw=resid2_raw,
           trend=trend2,
           label="Iron + cellulose prediction",
           color=GREEN,
           title="With hydrous material (cellulose) — excellent fit",
           chi2_label=f"Reduced χ² = {chi2_comb:.2f}",
           annotation=tape_note)

fig.suptitle(
    "Detecting ~100 µm Painter's Tape on a Pepper Can  |  Neutron TOF Transmission\n"
    "The full TOF spectrum reveals the hydrogen signature — invisible to white-beam radiography",
    fontsize=11.5, fontweight="bold",
)

outfile = "tape_detection_comparison.png"
plt.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
print(f"\nSaved: {outfile}")
