import { P } from './params.js';
import * as cardiac from './model.js';
import { initialState as cardiacInit, IDX as cardiacIDX } from './state.js';
import { rk4 } from './integrate.js';

/** Default model. Pass `model:` to runPopulation to use the HSPC one instead. */
const CARDIAC = { rhs: cardiac.rhs, initialState: cardiacInit, IDX: cardiacIDX };

// ---------------------------------------------------------------------------
// A dish is not one cell. The TBXT and GATA4 switches in this model are
// bistable, so a single deterministic cell gives an all-or-none answer, while
// flow cytometry gives you a percentage. The graded, bell-shaped dose-response
// people actually measure IS the heterogeneity: cells differ in how much drug
// they take up, how much destruction-complex capacity they carry, and how
// deeply pluripotent they were at t = 0, so they cross the switch threshold at
// different doses.
//
// Predicted consequence: the width of a line's CHIR window is set by its
// cell-to-cell variance, not by its mean sensitivity. Two lines with identical
// median IC50 but different CV will have very different "optimal dose" ranges.
// ---------------------------------------------------------------------------

/** Deterministic PRNG so every run is reproducible. */
export function rng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6D2B79F5) >>> 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Box-Muller, then exponentiate: lognormal with the given median and CV. */
function lognormal(rand, median, cv) {
  const sigma = Math.sqrt(Math.log(1 + cv * cv));
  const u1 = Math.max(rand(), 1e-12), u2 = rand();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return median * Math.exp(sigma * z);
}

/** Which parameters vary between cells, and by how much. */
export const DEFAULT_CV = {
  'chir.Kp':        0.25,  // drug uptake / efflux
  'wnt.k_phos':     0.30,  // destruction-complex capacity (Axin/APC/GSK3 dosage)
  'wnt.k_syn_b':    0.15,  // beta-catenin synthesis
  'wnt.k_sec_wnt':  0.35,  // local paracrine density - the biggest field effect
  'grn.k_tbxt':     0.20   // competence of the streak switch
};

const clone = (o) => JSON.parse(JSON.stringify(o));
const setPath = (o, path, v) => {
  const k = path.split('.');
  o[k[0]][k[1]] = v;
  if (k[1].startsWith('k_') || k[1].startsWith('d_')) {
    // GRN TFs keep k == d so they stay bounded on [0,1]
    const twin = (k[1].startsWith('k_') ? 'd_' : 'k_') + k[1].slice(2);
    if (k[0] === 'grn' && twin in o.grn) o.grn[twin] = v;
  }
};

export function makeCell(rand, cv = DEFAULT_CV) {
  const par = clone(P);
  for (const [path, c] of Object.entries(cv)) {
    const k = path.split('.');
    setPath(par, path, lognormal(rand, P[k[0]][k[1]], c));
  }
  return par;
}

/**
 * Run n independent cells through one protocol.
 * Returns the per-cell endpoint plus dish-level summaries.
 */
export function runPopulation(u, {
  n = 200, tEnd = 240, h = 0.01, seed = 1234, cv = DEFAULT_CV,
  markers = ['TNNT2', 'NKX25', 'GATA4', 'TBXT', 'SOX17', 'OCT4'],
  posThreshold = 0.5, model = CARDIAC
} = {}) {
  const { rhs, initialState, IDX } = model;
  const rand = rng(seed);
  const cells = [];
  for (let i = 0; i < n; i++) {
    const par = makeCell(rand, cv);
    const tr = rk4(rhs, initialState(), { tEnd, h, sample: tEnd, u, par });
    const y = tr[tr.length - 1].y;
    const rec = { V: y[IDX.V] };
    for (const mk of markers) rec[mk] = y[IDX[mk]];
    cells.push(rec);
  }

  const summary = { n, viability: mean(cells.map((c) => c.V)) };
  const totalMass = cells.reduce((a, c) => a + c.V, 0);
  for (const mk of markers) {
    // %positive is weighted by survival: dead cells are not in the flow gate
    summary[`pct_${mk}`] = 100 * cells.reduce(
      (a, c) => a + (c[mk] > posThreshold ? c.V : 0), 0) / Math.max(totalMass, 1e-9);
    summary[`mean_${mk}`] = mean(cells.map((c) => c[mk]));
  }
  // Yield per INPUT cell, not per surviving cell: what you actually harvest.
  // Only defined when the marker was requested, so this is skipped for models
  // that do not have it.
  for (const mk of markers) {
    summary[`yield_${mk}`] = cells.reduce(
      (a, c) => a + (c[mk] > posThreshold ? c.V : 0), 0) / n * 100;
  }
  return { cells, summary };
}

const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
