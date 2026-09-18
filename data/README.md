# Data format

One long-format CSV. One row per measurement. Put your own file here and point
`experiments/08_fit.js` at it.

```csv
kind,protocol,readout,time_h,value,sd,ref
```

| column | meaning |
|---|---|
| `kind` | `topflash`, `qpcr`, or `flow` |
| `protocol` | a name defined in the `PROTOCOLS` map of the fit script |
| `readout` | a model species — `TCF_B`, `TBXT`, `NKX25`, `TNNT2`, `HOXB4`, `CD34`, … |
| `time_h` | hours from the start of induction |
| `value` | see below — **fold-change** for qpcr/topflash, **percent** for flow |
| `sd` | optional. Fractional CV for fold-changes (0.35 = ±35 %), percentage points for flow |
| `ref` | optional. Another protocol name to take the fold-change against |

## What `value` has to be

**This is the part that goes wrong.** The three assays are compared differently.

- **`topflash`** — fold-change of the reporter relative to a reference. If `ref`
  is set, the model divides by that protocol at the same timepoint; if it is
  blank, it divides by the same protocol at t = 0. Do not put raw luminescence
  here; put the ratio you would plot.

- **`qpcr`** — the same, a fold-change (2^-ΔΔCt), not a Ct and not a raw
  normalised value. The model's transcription factors are on a 0–1 scale with a
  0.02 leak floor, so only ratios are comparable between model and bench.

- **`flow`** — percent positive, 0–100, of the **surviving** gate. This is the
  only readout scored against the population model rather than a single cell,
  so it is also the slowest to evaluate: put fewer of these in than qPCR points.

## Which parameters to fit with which data

Fit them in this order. Fitting everything at once against everything at once
produces a good SSE and meaningless parameters.

| Data you have | Fit | Leaves alone |
|---|---|---|
| TOPflash dose–response | `wnt.k_phos`, `wnt.k_syn_ax2`, `wnt.k_tcf_off`, `chir.Kp` | the whole fate layer |
| qPCR timecourse of TBXT/MESP1 | `grn.K_tcf_tbxt`, `grn.d_tbxt`, `grn.K_str` | the Wnt layer, once it is pinned |
| Flow % TNNT2+ vs dose | `population` CVs, `chir.IC50_tox` | everything already fitted |
| Cord blood qPCR (HOXB4, MYC) | `hsc.K_w_hoxb4`, `hsc.K_w_myc` | — |
| Cord blood CD34+ fold expansion | `hsc.k_prolif`, `hsc.K_p57` | — |

The two HSPC thresholds `K_w_hoxb4` and `K_w_myc` are the ones worth measuring
first: **their separation is what makes the cord blood dose-response
non-monotonic**, and nothing else in the model sets where the optimum sits.

## Identifiability

`fit()` constrains each parameter to within `bounds`× of its current value and
flags any that finish **at a bound** — that means your data does not constrain
it and the optimiser just pushed it as far as it was allowed. Treat a flagged
parameter as un-fitted, not as fitted-to-an-extreme.

Run `npm run fitcheck` first. It fits synthetic data generated from known
parameters and reports how well it recovers them. If recovery is poor on
synthetic data it will be worse on yours, and the fix is more informative
timepoints, not more iterations.

## Two things that will bite you

**Flow rows are ~60× the cost of every other row.** A `flow` row runs the whole
population ensemble, not one cell, so a handful of them turns a one-minute fit
into an hour. `09_fit_mydata.js` excludes them by default. That is a modelling
decision as much as a speed one: qPCR and TOPflash pin the signalling and fate
layers, while flow percentages mostly constrain the population CVs in
`src/population.js`. Fit them jointly and a CV will happily absorb error that
belongs to a rate constant. Fit in two passes instead.

**Which model a row belongs to is decided by its `protocol`, not its readout.**
`TCF_B` exists in both the cardiac and cord blood models, so filtering by
readout would quietly simulate an hPSC TOPflash experiment through the cord
blood model, which starts from completely different conditions. Declare each
protocol's cell type in the `PROTO_MODEL` map.

## Running a fit

```bash
node experiments/09_fit_mydata.js --data data/mine.csv --model cardiac
node experiments/09_fit_mydata.js --data data/mine.csv --model hsc
```

Add any protocol your rows refer to to the `PROTOCOLS` map at the top of that
file first — an unknown protocol name is a hard error rather than a silently
skipped row.
