"""
Python wrapper around the JavaScript model.

The model itself stays in JavaScript. This module only starts it and turns its
output into DataFrames, so there is exactly ONE implementation of the biology
and nothing to keep in sync. Every fix made to src/ is automatically reflected
here; a Python port would have needed each one done twice.

    import stemcell as sc

    df = sc.run(model="hsc", chir=1, chir_end=168, days=7)   # trajectory
    s  = sc.population(model="cardiac", chir=6, cells=200)   # flow-like summary
    sw = sc.sweep([0, 1, 2, 4, 6, 12], model="hsc", days=7, chir_end=168)

Requires Node on PATH. Nothing else.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

# The repo root, i.e. the directory holding run.js — this file lives in notebooks/
ROOT = Path(__file__).resolve().parent.parent
RUN_JS = ROOT / "run.js"


class ModelError(RuntimeError):
    """The model process failed. Carries its stderr, which is usually the answer."""


def _node() -> str:
    exe = shutil.which("node")
    if exe is None:
        raise ModelError(
            "node was not found on PATH. The model is JavaScript; install Node "
            "from https://nodejs.org and reopen the notebook so it picks up PATH."
        )
    return exe


def _invoke(**kwargs) -> dict:
    """Call run.js --json and parse stdout."""
    if not RUN_JS.exists():
        raise ModelError(f"cannot find {RUN_JS}. Is stemcell.py still inside notebooks/?")

    argv = [_node(), str(RUN_JS), "--json"]
    for key, value in kwargs.items():
        if value is None:
            continue
        argv += [f"--{key.replace('_', '-')}", str(value)]

    proc = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ModelError(
            f"model exited {proc.returncode}\n"
            f"command: {' '.join(argv)}\n"
            f"stderr:\n{proc.stderr.strip()}"
        )
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        # Almost always means something printed alongside the JSON.
        raise ModelError(
            f"could not parse model output as JSON: {exc}\n"
            f"first 400 chars:\n{proc.stdout[:400]}"
        ) from exc


def run(model: str = "cardiac", chir: float = 6, chir_start: float = 0,
        chir_end: float = 24, days: float = 10, iwp: str | None = None,
        activin: float | None = None, bmp4: float | None = None,
        wnt3a: float | None = None, cytokines: float | None = None) -> pd.DataFrame:
    """
    One deterministic cell. Returns a DataFrame indexed by time in hours, with
    one column per species plus GSK3_active.

    For the HSPC model the engraftment read-outs are attached as ``df.attrs``:
    fold_expansion, stem_frequency, homing, units, engraftment.
    """
    out = _invoke(model=model, chir=chir, chir_start=chir_start, chir_end=chir_end,
                  days=days, iwp=iwp, activin=activin, bmp4=bmp4, wnt3a=wnt3a,
                  cytokines=cytokines)
    df = pd.DataFrame(out["species"])
    df.insert(0, "GSK3_active", out["GSK3_active"])
    df.index = pd.Index(out["t_h"], name="t_h")
    df.attrs["protocol"] = out["protocol"]
    if "engraftment" in out:
        df.attrs.update(out["engraftment"])
    return df


def population(model: str = "cardiac", cells: int = 200, **kwargs) -> pd.Series:
    """
    A dish rather than a cell: an ensemble with realistic cell-to-cell spread.

    Returns the summary a flow cytometer would give you — pct_<marker> is the
    percentage of the SURVIVING gate above threshold, and yield_<marker> is per
    input cell. Use this, not run(), whenever you want to compare against flow.
    """
    out = _invoke(model=model, cells=cells, **kwargs)
    return pd.Series(out["summary"], name=f"{model}/{cells}cells")


def sweep(doses, model: str = "cardiac", cells: int = 1, **kwargs) -> pd.DataFrame:
    """
    Run a CHIR dose series. One row per dose.

    cells=1 gives deterministic endpoint levels; cells>1 gives population
    percentages, which is what you want for anything compared to flow data.
    """
    rows = []
    for dose in doses:
        if cells == 1:
            df = run(model=model, chir=dose, **kwargs)
            row = df.iloc[-1].to_dict()
            row.update({k: v for k, v in df.attrs.items() if k != "protocol"})
        else:
            row = population(model=model, chir=dose, cells=cells, **kwargs).to_dict()
        row["chir_uM"] = dose
        rows.append(row)
    return pd.DataFrame(rows).set_index("chir_uM")


def experiment(name: str) -> str:
    """
    Run one of the scripts in experiments/ and return its printed output.
    Use for the ones with no JSON mode, e.g. experiment("04_sensitivity").
    """
    path = ROOT / "experiments" / (name if name.endswith(".js") else name + ".js")
    if not path.exists():
        available = sorted(p.stem for p in (ROOT / "experiments").glob("*.js"))
        raise ModelError(f"no experiment {name!r}. Available: {', '.join(available)}")
    proc = subprocess.run([_node(), str(path)], cwd=ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ModelError(f"{name} exited {proc.returncode}\n{proc.stderr.strip()}")
    return proc.stdout


def load_csv(name: str) -> pd.DataFrame:
    """
    Read a CSV an experiment wrote into out/. Names are as written there, e.g.
    load_csv("06_cordblood_dose"). Run the experiment first if it is missing.
    """
    path = ROOT / "out" / (name if name.endswith(".csv") else name + ".csv")
    if not path.exists():
        existing = sorted(p.stem for p in (ROOT / "out").glob("*.csv"))
        raise ModelError(
            f"{path.name} does not exist. Run the experiment that writes it first.\n"
            f"Present in out/: {', '.join(existing) if existing else '(nothing yet)'}"
        )
    return pd.read_csv(path)
