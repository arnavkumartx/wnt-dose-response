"""
The model's headline figures, as importable functions and as a script.

In a notebook — every function returns a matplotlib Figure and also writes a
PNG, so the plot appears inline and the file lands in out/figures/:

    import make_figures as mf

    mf.fig_cordblood_dose()          # one figure
    figs = mf.all_figures()          # all of them, all drawn inline
    mf.fig_cardiac_bell(cells=200)   # same figure, bigger ensemble

From a shell — same code, headless, straight to files:

    python notebooks/make_figures.py --dpi 300

Every function takes dpi and nothing else it does not need, so they are
interchangeable. Shared simulation work is cached, so calling all of them costs
barely more than calling the slowest one.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

# Force the headless backend only when this file is run as a script. Switching
# it on import would stop inline plots rendering for the rest of a notebook
# session, with no error raised.
_AS_SCRIPT = __name__ == "__main__"
if _AS_SCRIPT:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

import stemcell as sc

OUT = sc.ROOT / "out" / "figures"

plt.rcParams.update({
    "figure.figsize": (8, 4.8),
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 10,
})

_CACHE: dict = {}


def cordblood_sweep() -> pd.DataFrame:
    """The 7-day cord blood dose series. Cached: two figures share it."""
    if "cb" not in _CACHE:
        doses = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12]
        _CACHE["cb"] = sc.sweep(doses, model="hsc", chir_end=168, days=7)
    return _CACHE["cb"]


def _save(fig, name: str, dpi: int):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    # As a script, close to free memory. In a notebook, leave it open so the
    # call also draws the plot instead of appearing to do nothing.
    if _AS_SCRIPT:
        plt.close(fig)
    print(f"  {path.relative_to(sc.ROOT)}")
    return fig


def fig_cordblood_dose(dpi: int = 200):
    """Cord blood: fewer cells, better cells. The headline result."""
    cb = cordblood_sweep()
    fig, ax1 = plt.subplots()
    ax1.plot(cb.index, cb["fold_expansion"], "o-", color="tab:blue")
    ax1.set_xlabel("CHIR99021 (µM)")
    ax1.set_ylabel("fold expansion", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.plot(cb.index, cb["engraftment"], "s-", color="tab:red")
    ax2.set_ylabel("engraftment per input cell", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")

    best = cb["engraftment"].idxmax()
    ax2.axvline(best, ls="--", color="grey", lw=1)
    ax2.annotate(f"optimum {best} µM", xy=(best, cb["engraftment"].max()),
                 xytext=(best + 1.0, cb["engraftment"].max() * 0.97), color="grey")
    ax1.set_title("Cord blood CD34+, 7 days: fewer cells, better cells")
    fig.text(0.5, -0.06, "High-dose CHIR is worse than no CHIR at all.",
             ha="center", fontsize=9, style="italic", color="0.35")
    return _save(fig, "01_cordblood_dose", dpi)


def fig_threshold_gap(dpi: int = 200):
    """Why that curve is a hill: two Wnt targets with different thresholds."""
    cb = cordblood_sweep()
    fig, ax = plt.subplots()
    ax.plot(cb.index, cb["HOXB4"], "o-", label="HOXB4  (self-renewal)")
    ax.plot(cb.index, cb["MYC"], "s-", label="MYC  (niche exit)")
    ax.plot(cb.index, cb["STEM"], "^-", color="black", label="stem frequency")
    ax.set_xlabel("CHIR99021 (µM)")
    ax.set_ylabel("level (0–1)")
    ax.legend()
    ax.set_title("The gap between these two thresholds sets the optimum")
    fig.text(0.5, -0.06,
             "Measure these two by qPCR first — nothing else sets where the peak sits.",
             ha="center", fontsize=9, style="italic", color="0.35")
    return _save(fig, "02_threshold_gap", dpi)


def fig_cardiac_biphasic(dpi: int = 200):
    """Wnt has to go on, then off."""
    arms = {
        "GiWi: CHIR d0–1, IWP d3–5": dict(chir=6, chir_end=24, iwp="72:120"),
        "CHIR d0–1 only": dict(chir=6, chir_end=24, iwp="none"),
        "CHIR continuous": dict(chir=6, chir_end=240, iwp="none"),
        "no CHIR": dict(chir=0, chir_end=24, iwp="none"),
    }
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.8))
    for label, kw in arms.items():
        d = sc.run(model="cardiac", days=10, **kw)
        axL.plot(d.index / 24, d["TCF_B"], label=label)
        axR.plot(d.index / 24, d["TNNT2"], label=label)
    axL.set_xlabel("day"); axL.set_ylabel("TCF·β-catenin (nM)")
    axL.set_title("Wnt pathway output")
    axR.set_xlabel("day"); axR.set_ylabel("TNNT2")
    axR.set_title("cardiomyocyte marker")
    axR.legend(fontsize=8)
    fig.suptitle("Continuous CHIR makes excellent cardiac mesoderm and zero cardiomyocytes",
                 y=1.02)
    return _save(fig, "03_cardiac_biphasic", dpi)


def fig_cardiac_bell(dpi: int = 200, cells: int = 80):
    """Rising limb is efficiency, falling limb is toxicity."""
    key = f"bell{cells}"
    if key not in _CACHE:
        doses = [0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20]
        rows = [sc.population(model="cardiac", chir=d, chir_end=24, iwp="72:120",
                              days=10, cells=cells) for d in doses]
        pop = pd.DataFrame(rows, index=doses)
        pop.index.name = "chir_uM"
        _CACHE[key] = pop
    pop = _CACHE[key]

    fig, ax = plt.subplots()
    ax.plot(pop.index, pop["pct_TNNT2"], "o-", label="% TNNT2+ of survivors")
    ax.plot(pop.index, pop["yield_TNNT2"], "s-", label="% per input cell")
    ax.plot(pop.index, pop["viability"] * 100, "^-", color="grey", label="viability %")
    ax.set_xlabel("CHIR99021 (µM)"); ax.set_ylabel("percent")
    ax.legend()
    ax.set_title(f"Cardiac dose–response, {cells} simulated cells per point")
    fig.text(0.5, -0.06, "The graded curve IS the cell-to-cell heterogeneity.",
             ha="center", fontsize=9, style="italic", color="0.35")
    return _save(fig, "04_cardiac_bell", dpi)


def fig_gsk3_occupancy(dpi: int = 200):
    """Why the dose is µM and not nM. Needs `npm run gsk3` to have written its CSV."""
    df = sc.load_csv("01_gsk3_occupancy")
    fig, ax = plt.subplots()
    ax.plot(df["chir_medium_uM"], df["pct_inhibited"], color="tab:purple")
    ax.axhline(50, ls=":", color="grey", lw=1)
    ic50 = df.loc[(df["pct_inhibited"] - 50).abs().idxmin(), "chir_medium_uM"]
    ax.axvline(ic50, ls="--", color="grey", lw=1)
    ax.annotate(f"medium IC50 ≈ {ic50:.1f} µM", xy=(ic50, 50),
                xytext=(ic50 + 1.5, 38), color="grey")
    ax.axvspan(3, 12, alpha=0.10, color="tab:green")
    ax.text(7.5, 12, "dose range used\nin protocols", ha="center", fontsize=9, color="0.35")
    ax.set_xlabel("CHIR99021 in the medium (µM)")
    ax.set_ylabel("% GSK3 inhibited")
    ax.set_title("A 9.9 nM inhibitor needs µM doses: ATP competition costs 201×")
    return _save(fig, "05_gsk3_occupancy", dpi)


#: name -> function, in the order they tell the story.
FIGURES = {
    "cordblood_dose": fig_cordblood_dose,
    "threshold_gap": fig_threshold_gap,
    "cardiac_biphasic": fig_cardiac_biphasic,
    "cardiac_bell": fig_cardiac_bell,
    "gsk3_occupancy": fig_gsk3_occupancy,
}


def all_figures(dpi: int = 200, cells: int | None = None) -> dict:
    """
    Every figure. Returns {name: Figure}.

    In a notebook each one is drawn inline as well as written to out/figures/.
    A figure whose input is missing is skipped with a note rather than killing
    the whole run.
    """
    made = {}
    for name, func in FIGURES.items():
        kw = {"dpi": dpi}
        if name == "cardiac_bell" and cells is not None:
            kw["cells"] = cells
        try:
            made[name] = func(**kw)
        except sc.ModelError as exc:
            print(f"  skipped {name}: {exc.args[0].splitlines()[0]}")
    return made


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dpi", type=int, default=200, help="200 for screen, 300 for print")
    ap.add_argument("--cells", type=int, default=80, help="ensemble size for the bell curve")
    args = ap.parse_args()

    print(f"writing figures to {OUT.relative_to(sc.ROOT)}/ at {args.dpi} dpi\n")
    all_figures(dpi=args.dpi, cells=args.cells)
    print("\ndone")


if __name__ == "__main__":
    main()
