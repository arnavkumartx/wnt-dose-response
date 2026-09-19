# Contributing

## Running it

Node 18+ is the only requirement for the model itself.

```bash
npm run fast     # the quick suite, ~1 min
npm run all      # everything, ~8 min
```

Python 3.12 with `jupyterlab pandas matplotlib` is needed only for the notebook
and figure generator.

```bash
npm run jupyter
```

## Before changing anything

Run `npm run calibrate` first. It checks numerical convergence by step halving
and prints the undifferentiated steady state against its targets. If those are
off, nothing downstream means anything.

After a change, `npm run biphasic` is the fastest regression test: the GiWi arm
should give TNNT2 ≈ 0.75 and the continuous-CHIR arm ≈ 0.02. Those two numbers
catch most breakages.

## Claims live in scripts, not comments

Every quantitative statement in the README is produced by something in
`experiments/`. If a change moves a number, re-run the relevant script and
update the README rather than leaving the prose stale.

Where a script prints a conclusion, it should derive it from the data it just
computed instead of hard-coding it. `07_iwr_vs_iwp.js` is the pattern: it scans
the dose axis, detects whether a rescue window exists, and prints whichever
conclusion is true.

## Parameters

Every parameter in `src/params.js` carries a provenance tag:

| tag | meaning |
|---|---|
| `[M]` | measured — from a published measurement |
| `[D]` | derived — computed from measured quantities |
| `[F]` | fitted — tuned to reproduce a published qualitative result |
| `[A]` | assumed — a placeholder |

Keep the tag accurate when changing a value. `npm run sens` ranks parameters by
normalised sensitivity; 66 of 112 move neither endpoint by more than 2 % per
10 % change, so check there before spending effort on one.

Transcription factors use `k == d` so each stays bounded on [0,1] and reads as a
fraction of maximal expression. The `k_*` values are generated from the `d_*`
ones at load time — add a `d_` entry, not both.

## Fitting

`npm run fitcheck` fits synthetic data generated from known parameters and
reports recovery. Run it before fitting real data; if recovery is poor on
synthetic data it will be worse on real data, and the fix is more informative
timepoints rather than more iterations.

See [data/README.md](data/README.md) for the CSV format and which parameters to
fit against which assay.

## Figures

`docs/figures/` holds the plots rendered in the README. Regenerate with:

```bash
python notebooks/make_figures.py --dpi 160
cp out/figures/*.png docs/figures/
```

Only commit updated figures when the underlying numbers have actually changed —
PNG diffs are noise otherwise.
