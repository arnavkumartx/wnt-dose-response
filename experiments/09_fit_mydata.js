import * as cardiac from '../src/model.js';
import * as cardiacState from '../src/state.js';
import * as hspc from '../src/hsc.js';
import { loadObservations, fit, writeOverrides, estimateCost } from '../src/fit.js';
import { table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// Fit the model to a CSV of your own measurements.
//
//   node experiments/09_fit_mydata.js --data data/mine.csv --model cardiac
//   node experiments/09_fit_mydata.js --data data/mine.csv --model hsc
//
// Run 08_fit_check.js first. If the fitter cannot recover known parameters from
// synthetic data it will not recover them from yours.
// ---------------------------------------------------------------------------

const argv = process.argv.slice(2);
const arg = (n, d) => {
  const i = argv.indexOf('--' + n);
  return i >= 0 && argv[i + 1] && !argv[i + 1].startsWith('--') ? argv[i + 1] : d;
};

const dataFile = arg('data', 'data/template.csv');
const which = arg('model', 'cardiac');

// Every protocol a data row may refer to. Add yours here.
const PROTOCOLS = {
  none:   [],
  chir2:  [{ from: 0, to: 24, set: { chir: 2 } }],
  chir4:  [{ from: 0, to: 24, set: { chir: 4 } }],
  chir6:  [{ from: 0, to: 24, set: { chir: 6 } }],
  chir12: [{ from: 0, to: 24, set: { chir: 12 } }],
  giwi4:  [{ from: 0, to: 24, set: { chir: 4 } },  { from: 72, to: 120, set: { iwp: 0.95 } }],
  giwi6:  [{ from: 0, to: 24, set: { chir: 6 } },  { from: 72, to: 120, set: { iwp: 0.95 } }],
  giwi12: [{ from: 0, to: 24, set: { chir: 12 } }, { from: 72, to: 120, set: { iwp: 0.95 } }],
  // cord blood: CHIR held for the whole 7-day culture
  cb_none:  [{ from: 0, to: 168, set: { cytokines: 1.0 } }],
  cb_chir1: [{ from: 0, to: 168, set: { chir: 1, cytokines: 1.0 } }],
  cb_chir6: [{ from: 0, to: 168, set: { chir: 6, cytokines: 1.0 } }]
};

// Which cell type each protocol was run in. This is what decides where a row
// belongs, NOT the readout: TCF_B exists in both models, so filtering by
// readout would quietly simulate an hPSC TOPflash experiment through the cord
// blood model, which has different starting conditions entirely.
const PROTO_MODEL = Object.fromEntries(
  Object.keys(PROTOCOLS).map((k) => [k, k.startsWith('cb_') ? 'hsc' : 'cardiac']));

const MODELS = {
  cardiac: {
    model: { rhs: cardiac.rhs, initialState: cardiacState.initialState, IDX: cardiacState.IDX },
    species: new Set(Object.keys(cardiacState.IDX)),
    // Wnt module first, then the fate thresholds it feeds.
    params: ['wnt.k_phos', 'wnt.k_syn_ax2', 'chir.Kp', 'grn.K_tcf_tbxt', 'grn.d_tbxt']
  },
  hsc: {
    model: { rhs: hspc.rhs, initialState: hspc.initialState, IDX: hspc.IDX },
    species: new Set(Object.keys(hspc.IDX)),
    // The separation of these two thresholds is what makes the cord blood
    // dose-response non-monotonic, so they are the ones worth pinning.
    params: ['hsc.K_w_hoxb4', 'hsc.K_w_myc', 'hsc.k_prolif', 'chir.Kp']
  }
};

const cfg = MODELS[which];
if (!cfg) { console.error(`unknown --model "${which}"; expected cardiac or hsc`); process.exit(1); }

const all = loadObservations(dataFile);
for (const o of all) {
  if (!(o.protocol in PROTOCOLS)) {
    console.error(`row refers to unknown protocol "${o.protocol}". `
      + `Add it to the PROTOCOLS map at the top of this file.`);
    process.exit(1);
  }
  if (!cfg.species.has(o.readout) && PROTO_MODEL[o.protocol] === which) {
    console.error(`readout "${o.readout}" is not a species of the ${which} model`);
    process.exit(1);
  }
}
const obs = all.filter((o) => PROTO_MODEL[o.protocol] === which
  && (!o.ref || PROTO_MODEL[o.ref] === which));
const skipped = all.length - obs.length;

console.log(`=== 09  Fit to ${dataFile} (${which} model) ===\n`);
console.log(`${obs.length} rows for this model`
  + (skipped ? `, ${skipped} skipped (run in the other cell type)` : '') + '\n');
if (!obs.length) { console.error('nothing to fit'); process.exit(1); }

if (dataFile.endsWith('template.csv')) {
  console.log('NOTE: data/template.csv holds PLACEHOLDER numbers, not measurements.');
  console.log('      Fitting the model to made-up data tells you nothing. This run');
  console.log('      only demonstrates that the pipeline works end to end.\n');
}

// A flow row costs popN simulations instead of one, so a handful of them turns
// a one-minute fit into an hour-long one. They are excluded by default, and
// that is a modelling decision as much as a speed one: qPCR and TOPflash pin
// the signalling and fate layers, while flow percentages mostly constrain the
// population CVs, which are a separate set of parameters. Fitting them jointly
// lets a CV absorb error that belongs to a rate constant.
const popN = Number(arg('cells', 60));
const withFlow = argv.includes('--with-flow');
const flowRows = obs.filter((o) => o.kind === 'flow');
const used = withFlow ? obs : obs.filter((o) => o.kind !== 'flow');

if (flowRows.length && !withFlow) {
  const cost = estimateCost(flowRows, popN);
  console.log(`Excluding ${flowRows.length} flow row(s): they would cost`
    + ` ~${Math.round(cost.simHoursPerEval / 1000)}k sim-hours per evaluation`
    + ` at ${popN} cells.`);
  console.log('Fit these in a second pass against the population CVs only');
  console.log('(see src/population.js DEFAULT_CV). Force with --with-flow.\n');
}
if (!used.length) { console.error('no rows left to fit'); process.exit(1); }

const t0 = Date.now();
const res = fit(cfg.model, PROTOCOLS, used, cfg.params, { bounds: 5, restarts: 2, popN });
console.log(table(['parameter', 'before', 'after', 'fold', 'flag'],
  Object.entries(res.values).map(([k, v]) => [k, r2(v.before, 4), r2(v.after, 4),
    r2(v.fold, 2) + 'x', v.atBound ? 'AT BOUND - not identified' : ''])));

console.log(`\nSSE  ${r2(res.sse_before, 3)}  ->  ${r2(res.sse_after, 3)}`
  + `   (${res.n} points, ${((Date.now() - t0) / 1000).toFixed(1)} s)`);

console.log('\nWorst-fitting points:');
console.log(table(['readout', 'protocol', 't (h)', 'measured', 'predicted', 'resid'],
  res.residuals.sort((a, b) => Math.abs(b.r) - Math.abs(a.r)).slice(0, 6)
    .map((x) => [x.readout, x.protocol, x.time_h, r2(x.value), r2(x.pred), r2(x.r)])));

const atBound = Object.entries(res.values).filter(([, v]) => v.atBound).map(([k]) => k);
if (atBound.length) {
  console.log(`\n${atBound.length} parameter(s) finished at a bound: ${atBound.join(', ')}`);
  console.log('Your data does not constrain these. Treat them as un-fitted and add');
  console.log('timepoints or conditions where they act differently from each other.');
}

const out = `out/overrides_${which}.js`;
writeOverrides(out, res.values);
console.log(`\n-> ${out}   (import OVERRIDES and applyOverrides(P, OVERRIDES) to use)`);
