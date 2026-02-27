"""
Tape Detection Demo – clean version
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
from scipy.ndimage import uniform_filter1d

RNG = np.random.default_rng(42)

# ── materials ──────────────────────────────────────────────────────────────
nbragg.register_material("Fe_sg229_Iron-alpha_CrysExtn1.ncmat")
nbragg.register_material("Cellulose_C6O5H10.ncmat")

# ── iron reference fit on pure data (no background, no extinction) ─────────
_raw = pd.read_csv("iron_powder_transmission.csv", sep=r"\s+",
                   header=0, names=["wl", "trans", "err"])
wl  = _raw["wl"].values
err = _raw["err"].values

xs_iron  = nbragg.CrossSection(iron=nbragg.materials["Fe_sg229_Iron-alpha_CrysExtn1.ncmat"])
model_fe = nbragg.TransmissionModel(xs_iron, vary_background=False, vary_extinction=False)

data_pure = nbragg.Data.from_transmission(_raw.rename(columns={"wl": "wavelength"}))
res_ref   = model_fe.fit(data_pure, wlmin=0.5, wlmax=9.0)
params_ref = res_ref.params.valuesdict()
T_iron_ref = model_fe.transmission(wl, **params_ref)

# ── cellulose (painter's tape) ─────────────────────────────────────────────
cel_mat = NC.load("Cellulose_C6O5H10.ncmat")
n_cel   = cel_mat.info.factor_macroscopic_xs
cel_xs  = (cel_mat.scatter.xsect(wl=wl, direction=(0,0,1)) +
           cel_mat.absorption.xsect(wl=wl, direction=(0,0,1)))

TAPE_TRUE_CM = 0.010    # 100 µm painter's tape
T_tape_true  = np.exp(-cel_xs * n_cel * TAPE_TRUE_CM)

# ── synthetic "pepper can + tape" dataset ─────────────────────────────────
T_data   = T_iron_ref * T_tape_true + RNG.normal(0., err)
err_data = err.copy()

df_tape   = pd.DataFrame({"wavelength": wl, "trans": T_data, "err": err_data})
data_tape = nbragg.Data.from_transmission(df_tape)

WLMIN, WLMAX = 0.5, 9.0
mask    = (wl >= WLMIN) & (wl <= WLMAX)
wl_fit  = wl[mask];  T_fit = T_data[mask];  err_fit = err_data[mask]
cel_fit = cel_xs[mask]

# ── Fit 1: iron only ───────────────────────────────────────────────────────
res_iron   = model_fe.fit(data_tape, wlmin=WLMIN, wlmax=WLMAX)
chi2_iron  = res_iron.redchi
T_iron_fit = model_fe.transmission(wl, **res_iron.params.valuesdict())
resid1     = (T_data - T_iron_fit) / err_data
print(f"Iron-only        χ²/dof = {chi2_iron:.3f}")

# ── Fit 2: iron + cellulose (joint, thickness free) ────────────────────────
# Start from reference iron params so optimizer can cleanly find t_cel
params_joint = lmfit.Parameters()
for k, v in params_ref.items():
    pobj = res_ref.params[k]
    params_joint.add(k, value=v, min=pobj.min, max=pobj.max, vary=pobj.vary)
params_joint.add("t_cel", value=TAPE_TRUE_CM * 0.5,   # start at 50 µm
                 min=1e-5, max=0.5, vary=True)

def resid_joint(p):
    kws  = {k: p[k].value for k in p if k != "t_cel"}
    T_fe = model_fe.transmission(wl_fit, **kws)
    T_ce = np.exp(-cel_fit * n_cel * p["t_cel"].value)
    return (T_fit - T_fe * T_ce) / err_fit

res_joint  = lmfit.Minimizer(resid_joint, params_joint).minimize(method="leastsq")
t_cel_fit  = res_joint.params["t_cel"].value
chi2_joint = res_joint.redchi
iron_best  = {k: res_joint.params[k].value for k in res_joint.params if k != "t_cel"}
T_joint_fit = (model_fe.transmission(wl, **iron_best) *
               np.exp(-cel_xs * n_cel * t_cel_fit))
resid2 = (T_data - T_joint_fit) / err_data
print(f"Iron + cellulose χ²/dof = {chi2_joint:.3f}")
print(f"Fitted tape thickness   = {t_cel_fit*1e4:.0f} µm  (true: {TAPE_TRUE_CM*1e4:.0f} µm)")

# ── smoothed residual trend (200-pt window) ────────────────────────────────
SMOOTH = 200
trend1 = uniform_filter1d(resid1, SMOOTH, mode="nearest")
trend2 = uniform_filter1d(resid2, SMOOTH, mode="nearest")

# ── figure ─────────────────────────────────────────────────────────────────
RED   = "#c0392b"
GREEN = "#27ae60"
GREY  = "0.4"

fig = plt.figure(figsize=(12, 7))
fig.patch.set_facecolor("white")
gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1], hspace=0.06, wspace=0.28,
                       left=0.08, right=0.97, top=0.91, bottom=0.09)
ax_l  = fig.add_subplot(gs[0, 0])
ax_lr = fig.add_subplot(gs[1, 0], sharex=ax_l)
ax_r  = fig.add_subplot(gs[0, 1])
ax_rr = fig.add_subplot(gs[1, 1], sharex=ax_r)

step = max(1, len(wl) // 400)

def panel(ax_m, ax_r_, T_fit_, resid_, trend_, color, title, chi2, note=None):
    ax_m.errorbar(wl[::step], T_data[::step], yerr=err_data[::step],
                  fmt=".", color=GREY, ms=2.5, lw=0.4, alpha=0.5, zorder=1)
    ax_m.plot(wl, T_fit_, color=color, lw=2.0, zorder=3)
    ax_m.set_xlim(WLMIN, WLMAX);  ax_m.set_ylabel("Transmission", fontsize=10)
    ax_m.tick_params(labelbottom=False)
    ax_m.set_title(title, fontsize=11, fontweight="bold")
    ax_m.text(0.03, 0.06, f"χ²/dof = {chi2:.2f}",
              transform=ax_m.transAxes, fontsize=9,
              bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9))
    if note:
        ax_m.text(0.97, 0.06, note, transform=ax_m.transAxes,
                  ha="right", va="bottom", fontsize=9, color=color,
                  bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=color, alpha=0.9))

    ax_r_.axhline(0, color="0.35", lw=0.8, ls="--")
    ax_r_.fill_between(wl, -1, 1, color="#d5f5e3", alpha=0.6)
    ax_r_.plot(wl, resid_, ".", color=color, ms=1.2, alpha=0.22)
    ax_r_.plot(wl, trend_,  color=color, lw=2.0, alpha=0.9)
    ax_r_.set_ylabel("Residual (σ)", fontsize=9)
    ax_r_.set_xlabel("Wavelength (Å)", fontsize=10)
    ax_r_.set_xlim(WLMIN, WLMAX)
    ylim = max(3.5, np.percentile(np.abs(resid_[mask]), 99) * 0.65)
    ax_r_.set_ylim(-ylim, ylim)
    ax_r_.yaxis.set_major_locator(plt.MultipleLocator(2))

panel(ax_l, ax_lr, T_iron_fit, resid1, trend1, RED,
      "Without hydrous material", chi2_iron)

panel(ax_r, ax_rr, T_joint_fit, resid2, trend2, GREEN,
      "With hydrous material (cellulose)",  chi2_joint,
      note=f"tape thickness: {t_cel_fit*1e4:.0f} µm")

fig.suptitle("Detecting painter's tape (~100 µm cellulose) on a pepper can\n"
             "Neutron TOF transmission — sensitivity enabled by the full wavelength spectrum",
             fontsize=11, fontweight="bold")

outfile = "tape_detection_comparison.png"
plt.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
print(f"Saved {outfile}")
