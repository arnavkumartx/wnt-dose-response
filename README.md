# CHIR99021 → GSK3 → Wnt/β-catenin → stem cell fate

A mechanistic, runnable model of how a small molecule at the top of a signalling
cascade ends up choosing a cell fate. Zero dependencies, pure Node.

```bash
npm run fast     # calibration, GSK3 curve, biphasic test, cord blood, inhibitors  (~1 min)
npm run all      # adds the population dose-response, lineage map and sensitivity  (~8 min)
```

Two fate modules sit on one shared signalling core:

| module | file | what it is for |
|---|---|---|
| **Cardiac mesoderm** | `src/model.js` | hPSC → cardiomyocyte, the canonical CHIR/GiWi protocol |
| **Cord blood HSPC** | `src/hsc.js` | CD34+ ex vivo culture, self-renewal vs engraftment |

---

## Read this first if you are working on cord blood

The model says — and the literature agrees — that **in cord blood HSPC, GSK3
inhibition is a self-renewal and engraftment intervention, not a
differentiation one.** According to PubMed:

- [Holmes 2008, *Stem Cells*](https://doi.org/10.1634/stemcells.2007-0600) —
  GSK-3β inhibition in cord blood CD34+ raises β-catenin, c-myc and HoxB4;
  it **delays** expansion while **preserving** stem cell activity, and raises
  CXCR4 and stromal adherence.
- [Ko 2011, *Stem Cells*](https://doi.org/10.1002/stem.551) — lengthens cycle
  time, raises p57, lowers cyclin D1, and improves engraftment *per repopulating
  unit* rather than unit number.
- [Luis 2011, *Cell Stem Cell*](https://doi.org/10.1016/j.stem.2011.07.017) —
  canonical Wnt regulates haematopoiesis in a **dosage-dependent** way. Five
  graded Apc alleles give lineage-specific optima, which is why gain-of-function
  studies report HSC expansion *or* HSC depletion depending on how hard they
  pushed.
- [Shen 2014, *Stem Cells Dev*](https://doi.org/10.1089/scd.2014.0230) — the
  same inhibitor given *in vivo* after transplant **impairs** regeneration and
  causes lasting T-cell depletion.

So pushing CHIR harder to drive differentiation is the wrong direction for this
cell type. The model is built to show that rather than to argue it: run
`npm run cordblood`.

A regulatory note, since it motivated the cell choice: cord blood units and the
expanded product omidubicel are FDA-approved, but that approval does not extend
to CHIR99021, which is a research compound with no approved human use. GSK3
inhibitors for HSPC remain preclinical.

---

## The four layers

```
  CHIR99021 (medium, µM)
        │  uptake, Kp
        ▼
  GSK3α/β activity          Ki_app = Ki·(1 + [ATP]/Km_ATP)   ← ATP-competitive
        │
        ▼
  destruction complex   ◄── AXIN2 (negative feedback, a Wnt target gene)
  Axin·APC·GSK3·CK1α    ◄── Fzd/LRP6 ◄── WNT3A / autocrine WNT ⊣ IWP2
        │               ◄── Axin stabilised by IWR-1 (tankyrase)
        ▼
  β-catenin:  cytoplasmic ⇄ nuclear ⇄ TCF/LEF     ⇄ cadherin pool
        │                      ▲
        │                      └── ARM-groove competition: TCF vs ICAT vs APC
        ▼                          BCL9/Pygo binds ARM1 — required, not competing
  TCF·β-catenin  ── the pathway output
        │
        ├──► CARDIAC:  TBXT → MIXL1 → MESP1 → GATA4 → NKX2-5 → TNNT2
        │              EOMES → SOX17 (needs Nodal)
        │              ⊣ SOX2/NANOG/OCT4 ;  TCF·β-cat ⊣ NKX2-5
        │
        └──► HSPC:     HOXB4 (low threshold)  → self-renewal
                       MYC   (high threshold) → niche exit, differentiation
                       p57 ↑, cyclin D1 ↓, CXCR4 ↑
```

### Layer 0 — why 6 µM and not 6 nM

CHIR99021 is **ATP-competitive** with a ~10 nM biochemical Ki. Cells hold 1–5 mM
ATP and GSK3β has a 15 µM Km for it, so Cheng–Prusoff gives

```
Ki_app = Ki · (1 + [ATP]/Km_ATP) = 9.9 nM · (1 + 3000/15) ≈ 2.0 µM
```

a **201-fold** shift, ≈3.6 µM once divided by the cell:medium partition
coefficient. That is the entire explanation for the dose every protocol uses,
and it comes out of measured numbers with no fitting.

It also predicts that anything moving the cellular ATP pool moves the effective
dose: across 1–5 mM the medium IC50 slides 1.2 → 6.0 µM. Galactose medium,
hypoxia, oligomycin, or plating density should all shift your dose–response — a
candidate explanation for batch-to-batch variability that is testable in a day.

### Layer 1 — β-catenin

Lee et al. 2003 architecture: Axin is the limiting scaffold, the destruction
complex phosphorylates β-catenin with Michaelis–Menten kinetics, and AXIN2 is
itself a Wnt target, closing a strong negative feedback loop. Added on top:

- **Nuclear retention.** APC and Axin shuttle β-catenin *out* of the nucleus, so
  inhibiting the complex both spares β-catenin and retains it — two gains in
  series, which is why the TCF response is steeper than total β-catenin.
- **The cadherin sink.** In hPSC most β-catenin is at junctions. EMT at the
  primitive streak releases it: feedback from fate back onto signalling.
- **Explicit scaffold assembly** (`src/interfaces.js`), so a tankyrase inhibitor
  and a porcupine inhibitor are different perturbations rather than one knob.

### Layer 2 — fate

Every transcription factor obeys `dX/dt = k·drive − d·X` with `k == d`, so each
is bounded on [0,1] and reads as "fraction of maximal expression".

Load-bearing wiring, cardiac:

| Wiring | Why |
|---|---|
| `TCF·β-cat → TBXT`, gated by SOX2 | the pluripotency-exit switch |
| SOX2 *raises TBXT's threshold*, does not veto it | a hard veto makes the switch unflippable |
| `S(Activin) AND TCF·β-cat → EOMES → SOX17` | the mesoderm/endoderm fork |
| `TCF·β-cat ⊣ NKX2-5` | forces the protocol to be biphasic |
| GATA4 self-activation | carries cardiac identity across the Wnt-off gap |
| SOX17/FOXA2 self-activation | the same job on the endoderm side |
| NKX2-5 self-activation | commitment; returning Wnt no longer reverses it |

**Both lineages need a latch to survive the CHIR washout.** Remove GATA4's
self-activation and the model makes no cardiomyocytes on any schedule; remove
SOX17's and endoderm peaks on day 1 and collapses when the 24 h spike ends.
That symmetry is a prediction: the commitment step in each arm should be a
self-reinforcing transcriptional loop, not a downstream response to Wnt.

Load-bearing wiring, HSPC — **the non-monotonic dose-response is not asserted
anywhere.** It emerges because two direct Wnt targets have different thresholds:
HOXB4 (self-renewal) switches on *below* MYC (niche exit). Mild activation
raises self-renewal; strong activation overruns it.

### Layer 3 — the dish

A dish is not one cell. The TBXT and GATA4 switches are bistable, so one
deterministic cell gives an all-or-none answer while flow cytometry gives a
percentage. `population.js` runs an ensemble with log-normal spread in drug
uptake, destruction-complex capacity, β-catenin synthesis, local paracrine
density and streak competence. **The graded dose–response curve you measure is
the heterogeneity.**

---

## Results

### The biphasic Wnt requirement falls out (`npm run biphasic`)

Nothing was fitted to this. It follows from one repression term plus two latches.

| arm | peak TBXT | GATA4 d10 | NKX2-5 d10 | TNNT2 d10 |
|---|---|---|---|---|
| no CHIR | 0.05 | 0.03 | 0.02 | **0.02** |
| CHIR 6 µM d0–d1 only | 0.91 | 0.95 | 0.11 | **0.13** |
| CHIR 6 µM continuous | 0.96 | 0.98 | 0.02 | **0.02** |
| GiWi (CHIR d0–1, IWP d3–5) | 0.91 | 0.87 | 0.72 | **0.75** |
| IWP d3–5 only | 0.05 | 0.03 | 0.02 | **0.02** |

Continuous CHIR makes excellent cardiac *mesoderm* and zero cardiomyocytes.

### Cord blood: the optimum is interior and low (`npm run cordblood`)

7-day CD34+ culture with SCF/TPO/FLT3L:

| CHIR | fold expansion | stem frequency | engraftment |
|---|---|---|---|
| 0 µM | 48.0× | 0.119 | 2.74 |
| **1 µM** | **30.6×** | **0.262** | **4.96**  (1.81× control) |
| 6 µM | 5.5× | 0.202 | 0.71 |
| 12 µM | 1.4× | 0.199 | ~0 |

Fewer cells, better cells, and past the MYC threshold neither. **High-dose CHIR
is worse than no CHIR at all.** Lymphoid output falls monotonically with dose —
there is no lymphoid-friendly window anywhere.

### Bell-shaped cardiac dose–response (`npm run dose`)

Rising limb = differentiation efficiency, falling limb = viability. Optimum
6 µM at 64.4 % yield, >50 %-of-max window 3–16 µM, and the top is fairly flat
(8 µM gives 61.9 %). Scaling the cell-to-cell CV trades peak yield against
window width: **a tighter clone peaks higher but has a narrower usable dose
range.**

One honest wart: ~28 % of cells are still TBXT+ at day 10 at every dose above
3 µM. These are the cells with the highest drawn autocrine Wnt tone, which never
get a Wnt-low window and so never commit or extinguish the streak programme.
That is a prediction, not a fudge — **the non-cardiomyocyte fraction should be
TBXT+ residual mesoderm concentrated where local Wnt was highest** — but it is
higher than a real day-10 culture and is a place the model is probably wrong.

### Wnt and Nodal are two separate knobs (`npm run lineage`)

At 26 h (where TBXT and MESP1 peak) and day 7 for SOX17, across CHIR × Activin:

| | Activin 0.18 (autocrine floor) | Activin 1.0 (100 ng/mL) |
|---|---|---|
| **CHIR 0** | nothing leaves pluripotency | 79 % SOX17+ |
| **CHIR 9** | 98 % TBXT+, 97 % MESP1+ | 100 % SOX17+, MESP1 ~0 |

Wnt sets whether cells leave pluripotency at all; Nodal sets which side of the
streak they leave through. Neither alone picks a lineage, because EOMES is an
AND-gate on both.

Two cautions this experiment taught the hard way, and both apply to your bench
readouts as much as to the model:

- **MESP1 is transient.** It peaks at 24–26 h, is down to ~0.3 by 48 h, and is
  gone by day 3. Scored at day 7 it reads as background at every dose and looks
  like the cardiac arm simply does not work.
- **GATA4 does not discriminate these lineages.** GATA4/GATA6 are definitive
  endoderm TFs as well as cardiac mesoderm ones, and it goes high in both
  corners here. MESP1 is the cardiac-specific call.

Every panel reports mean level next to % positive, because MESP1 peaks at ~0.57
against a 0.5 gate — the percentage there swings on small parameter changes
while the mean does not, and reporting only the percentage would let a gating
artefact read as a biological effect.

### What the interface layer does and does not predict (`npm run inhibitors`)

With CHIR withdrawn, IWP2 and IWR-1 give the same fate. With CHIR still on,
neither rescues — a 2.35× scaffold increase cannot outrun direct GSK3
inhibition. Where they *do* separate is the depth of the Wnt drop during the
washout window (TCF_B 0.82 for IWR-1 vs 1.55 for IWP2). **Testable and cheap:
run GiWi with each inhibitor and read AXIN2 by qPCR at d4.**

ICAT titration shows the other half: ICAT and TCF are mutually exclusive on ARM
3–9, so ICAT drops the output ~60 % without changing β-catenin levels at all —
a western blot would look unchanged while the reporter falls.

---

## Fitting it to your data

```bash
npm run fitcheck     # recovers known parameters from synthetic data first
```

See [data/README.md](data/README.md) for the CSV format and — more importantly —
which parameters to fit against which assay. Fitting everything against
everything gives a good SSE and meaningless parameters.

`fit()` bounds each parameter within a multiplicative factor of its current
value and **flags any that finish at a bound**: that means your data does not
constrain it, not that it fitted to an extreme.

## Parameter provenance

Every parameter in `src/params.js` is tagged `[M]` measured, `[D]` derived,
`[F]` fitted to a published qualitative result, or `[A]` assumed. `npm run sens`
ranks all of them by normalised sensitivity, so you can see which guesses matter
before spending bench time.

The two thresholds worth measuring first are `hsc.K_w_hoxb4` and `hsc.K_w_myc`.
**Their separation is what makes the cord blood dose-response non-monotonic**,
and nothing else in the model sets where the optimum sits.

Two results from `npm run sens` that change how you should calibrate:

**66 of 112 parameters move neither endpoint by >2 % per 10 % change.** Leave
them at their assumed values. The whole top of the ranking is the Wnt module.

**The four CHIR pharmacology parameters are structurally unidentifiable
individually**, and invisible to a fate readout at all:

| parameter | S[TNNT2 d10] | S[peak TCF_B] |
|---|---|---|
| `chir.Ki_gsk3b` | 0.00 | −0.44 |
| `chir.ATP_cell` | 0.00 | −0.44 |
| `chir.Kp` | 0.00 | +0.44 |
| `chir.Km_ATP` | 0.00 | +0.43 |

They only ever appear as the combination `Ki·(1 + ATP/Km)/Kp`, so the equal and
opposite sensitivities are exact — only the apparent IC50 is identifiable, and
fitting them separately will return whatever the optimiser wandered into. The
zeros matter too: at 6 µM the system is saturated past the switch, so **TNNT2
tells you nothing about CHIR potency.** Use TOPflash, or a fate readout near the
threshold dose (2–3 µM), where the switch is still responsive.

## What this model is not

- **Not spatial.** Colony density and edge effects drive much of the patchiness
  in real differentiation. Every cell here sees the same medium.
- **Not structural.** The interface layer is binding thermodynamics, not atomic
  detail. It says nothing about how CHIR sits in the GSK3β ATP pocket.
- **Not niche-aware.** The HSPC module has no stroma, no hypoxia, no adhesion —
  CXCR4 is a number, not a marrow.
- **Not fitted to your data.** The `[F]` values reproduce published qualitative
  behaviour. Quantitative prediction for your line needs your numbers.

## Layout

```
src/params.js      every parameter, tagged by provenance
src/wnt-core.js    layers 0-1, shared by both fate modules
src/interfaces.js  destruction-complex assembly + ARM-groove competition
src/model.js       cardiac mesoderm fate layer
src/hsc.js         cord blood HSPC fate layer
src/population.js  heterogeneous ensemble -> flow-cytometry-like readouts
src/fit.js         Nelder-Mead calibration against your measurements
src/kinetics.js    Hill / Michaelis-Menten / noisy-OR helpers
src/integrate.js   RK4, step verified by halving, throws on non-finite states
src/protocol.js    medium schedules (GiWi, endoderm, cord blood)
experiments/       runnable studies, each writes CSV to out/
data/              where your measurements go
run.js             CLI, both models:
                     node run.js --chir 6 --chir-end 24 --days 10
                     node run.js --model hsc --chir 1 --chir-end 168 --days 7
```
