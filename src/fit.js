import fs from 'node:fs';
import { P } from './params.js';
import { rk4 } from './integrate.js';
import { protocol } from './protocol.js';
import { runPopulation } from './population.js';

// ---------------------------------------------------------------------------
// Fit the [F]- and [A]-tagged parameters to your own measurements.
//
// Three readout kinds, because the three assays measure different things and
// comparing them naively is the usual way a fit goes wrong:
//
//   topflash  fold-change of the reporter vs untreated  -> model TCF_B ratio
//   qpcr      fold-change vs a reference condition      -> model TF ratio
//   flow      % of cells above a gate                   -> POPULATION model
//
// Fold-changes are compared in log space: a 10x measurement and a 0.1x
// measurement are equally informative, and linear residuals would let the
// large numbers swamp the fit.
// ---------------------------------------------------------------------------

/** Parse the long-format CSV described in data/README.md. */
export function loadObservations(file) {
  const raw = fs.readFileSync(file, 'utf8').trim().split(/\r?\n/);
  // The header is the first line that is neither blank nor a comment - a data
  // file that opens with a comment block is the normal case, not an edge case.
  const headIdx = raw.findIndex((l) => l.trim() && !l.trimStart().startsWith('#'));
  if (headIdx < 0) throw new Error(`${file} has no header row`);
  const cols = raw[headIdx].split(',').map((s) => s.trim());
  for (const required of ['kind', 'protocol', 'readout', 'time_h', 'value']) {
    if (!cols.includes(required)) {
      throw new Error(`${file}: header is missing the "${required}" column `
        + `(found: ${cols.join(', ')})`);
    }
  }
  const lines = raw.slice(headIdx + 1);
  return lines.filter((l) => l.trim() && !l.trimStart().startsWith('#')).map((l) => {
    const v = l.split(',').map((s) => s.trim());
    const o = Object.fromEntries(cols.map((c, i) => [c, v[i]]));
    return {
      kind: o.kind, protocol: o.protocol, readout: o.readout,
      time_h: Number(o.time_h), value: Number(o.value),
      sd: o.sd ? Number(o.sd) : null,
      ref: o.ref || null            // reference condition for a fold-change
    };
  });
}

const getPath = (o, p) => p.split('.').reduce((a, k) => a[k], o);
const clone = (o) => JSON.parse(JSON.stringify(o));

function setPath(o, p, v) {
  const k = p.split('.');
  o[k[0]][k[1]] = v;
  // Keep every TF synthesis/degradation pair matched so it stays on [0,1].
  for (const grp of ['grn', 'hsc']) {
    if (k[0] !== grp) continue;
    const twin = k[1].startsWith('k_') ? 'd_' + k[1].slice(2)
               : k[1].startsWith('d_') ? 'k_' + k[1].slice(2) : null;
    if (twin && twin in o[grp]) o[grp][twin] = v;
  }
}

/**
 * Build the simulator the objective calls.
 *
 * @param model      { rhs, initialState, IDX }  cardiac (model.js) or HSPC (hsc.js)
 * @param protocols  { name: [window, ...] }     named medium schedules
 */
export function makeSimulator(model, protocols, { popN = 120, h = 0.02 } = {}) {
  const cache = new Map();
  return {
    reset: () => cache.clear(),

    /** Deterministic level of one readout at one time. */
    level(par, protoName, readout, timeH) {
      if (timeH <= 0) return model.initialState()[model.IDX[readout]];
      const key = protoName + '|' + timeH + '|det';
      let tr = cache.get(key);
      if (!tr) {
        const u = protocol(protocols[protoName]);
        tr = rk4(model.rhs, model.initialState(), { tEnd: timeH, h, sample: timeH, u, par });
        cache.set(key, tr);
      }
      return tr[tr.length - 1].y[model.IDX[readout]];
    },

    /** Population % positive: what a flow gate actually reports. */
    pct(par, protoName, readout, timeH) {
      const key = protoName + '|' + timeH + '|pop|' + readout;
      let s = cache.get(key);
      if (!s) {
        s = runPopulation(protocol(protocols[protoName]),
          { n: popN, tEnd: timeH, h, markers: [readout], seed: 4242 }).summary;
        cache.set(key, s);
      }
      return s['pct_' + readout];
    }
  };
}

/** Weighted residuals. Fold-changes in log space, percentages linearly. */
export function objective(par, obs, sim) {
  let sse = 0;
  const residuals = [];
  for (const o of obs) {
    if (o.kind === 'flow') {
      const pred = sim.pct(par, o.protocol, o.readout, o.time_h);
      const r = (pred - o.value) / (o.sd ?? 10);
      sse += r * r;
      residuals.push({ ...o, pred, r });
      continue;
    }
    const num = sim.level(par, o.protocol, o.readout, o.time_h);
    const den = o.ref ? sim.level(par, o.ref, o.readout, o.time_h)
                      : sim.level(par, o.protocol, o.readout, 0);
    // A marker that starts at zero (TBXT, NKX25, TNNT2 ...) has no meaningful
    // fold-change against t = 0. That is not a modelling subtlety - it is what
    // ddCt does too, which is why qPCR is normalised to a reference CONDITION.
    if (!(den > 1e-4)) {
      throw new Error(`fold-change for ${o.readout} has a ~zero denominator `
        + `(${o.ref ? 'ref protocol ' + o.ref : 't=0'}). Set "ref" on this row to an `
        + `undifferentiated control protocol instead.`);
    }
    const pred = num / den;
    const r = (Math.log(Math.max(pred, 1e-9)) - Math.log(Math.max(o.value, 1e-9)))
              / Math.log(1 + (o.sd ?? 0.35));
    sse += r * r;
    residuals.push({ ...o, pred, r });
  }
  return { sse, n: obs.length, residuals };
}

// --- Nelder-Mead -----------------------------------------------------------
// Derivative-free, no dependencies, and robust to the mildly discontinuous
// surface a bistable switch produces. Parameters are searched in log space, so
// each is explored multiplicatively and none can go negative.
export function nelderMead(f, x0, { maxIter = 400, tol = 1e-6, step = 0.25 } = {}) {
  const n = x0.length;
  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) {
    const p = x0.slice();
    p[i] += step;
    simplex.push(p);
  }
  let fv = simplex.map(f);

  const centroidExcluding = (ex) => {
    const c = new Array(n).fill(0);
    simplex.forEach((p, i) => {
      if (i !== ex) for (let j = 0; j < n; j++) c[j] += p[j] / n;
    });
    return c;
  };

  for (let iter = 0; iter < maxIter; iter++) {
    const order = fv.map((v, i) => i).sort((a, b) => fv[a] - fv[b]);
    simplex = order.map((i) => simplex[i]);
    fv = order.map((i) => fv[i]);
    if (Math.abs(fv[n] - fv[0]) < tol * (Math.abs(fv[0]) + tol)) break;

    const c = centroidExcluding(n);
    const xr = c.map((v, j) => v + (v - simplex[n][j]));            // reflect
    const fr = f(xr);
    if (fr < fv[0]) {
      const xe = c.map((v, j) => v + 2 * (v - simplex[n][j]));      // expand
      const fe = f(xe);
      if (fe < fr) { simplex[n] = xe; fv[n] = fe; } else { simplex[n] = xr; fv[n] = fr; }
    } else if (fr < fv[n - 1]) {
      simplex[n] = xr; fv[n] = fr;
    } else {
      const xc = c.map((v, j) => v + 0.5 * (simplex[n][j] - v));    // contract
      const fc = f(xc);
      if (fc < fv[n]) {
        simplex[n] = xc; fv[n] = fc;
      } else {                                                      // shrink
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((v, j) => simplex[0][j] + 0.5 * (v - simplex[0][j]));
          fv[i] = f(simplex[i]);
        }
      }
    }
  }
  const best = fv.map((v, i) => i).reduce((a, b) => (fv[a] < fv[b] ? a : b));
  return { x: simplex[best], f: fv[best] };
}

/**
 * Fit `paramPaths` to `obs`.
 *
 * `bounds` is a multiplicative factor around each current value: a parameter is
 * only allowed to move within [v/bounds, v*bounds]. If a fitted value lands on
 * a bound, that parameter is not identified by the data you supplied - the fit
 * reports it rather than hiding it.
 */
/**
 * Rough cost of one objective evaluation, in single-cell simulation-hours.
 * A `flow` row costs popN times what a deterministic row costs, which is the
 * difference between a fit that takes a minute and one that takes an hour.
 */
export function estimateCost(obs, popN = 60) {
  let hours = 0;
  for (const o of obs) {
    hours += o.time_h * (o.kind === 'flow' ? popN : 1);
    if (o.ref) hours += o.time_h;
  }
  return { simHoursPerEval: hours, flowRows: obs.filter((o) => o.kind === 'flow').length };
}

export function fit(model, protocols, obs, paramPaths, {
  bounds = 5, restarts = 3, popN = 60, seedNoise = 0.3
} = {}) {
  const sim = makeSimulator(model, protocols, { popN });
  const base = paramPaths.map((p) => getPath(P, p));
  const lo = base.map((v) => Math.log(v / bounds));
  const hi = base.map((v) => Math.log(v * bounds));

  // Validate once, loudly, with the starting parameters. Inside the optimiser a
  // failed run is caught and scored 1e12 (some parameter combinations genuinely
  // blow up), so a STRUCTURAL problem with the observations would otherwise be
  // swallowed and the fit would silently return its starting point.
  sim.reset();
  objective(P, obs, sim);

  const parFrom = (logx) => {
    const par = clone(P);
    logx.forEach((lv, i) => {
      setPath(par, paramPaths[i], Math.exp(Math.min(hi[i], Math.max(lo[i], lv))));
    });
    return par;
  };
  const evalAt = (logx) => {
    sim.reset();
    try { return objective(parFrom(logx), obs, sim).sse; }
    catch { return 1e12; }          // a non-finite run is simply rejected
  };

  const start0 = base.map(Math.log);
  let best = { f: Infinity, x: start0 };
  for (let r = 0; r < restarts; r++) {
    const s = r === 0 ? start0
      : start0.map((v, i) => v + (Math.random() * 2 - 1) * seedNoise * (hi[i] - lo[i]) / 2);
    const res = nelderMead(evalAt, s, { maxIter: 400 });
    if (res.f < best.f) best = res;
  }

  const fitted = parFrom(best.x);
  const values = {};
  best.x.forEach((lv, i) => {
    const v = Math.exp(Math.min(hi[i], Math.max(lo[i], lv)));
    const atBound = Math.abs(Math.log(v / base[i])) > 0.98 * Math.log(bounds);
    values[paramPaths[i]] = { before: base[i], after: v, fold: v / base[i], atBound };
  });

  sim.reset();
  const before = objective(P, obs, sim);
  sim.reset();
  const after = objective(fitted, obs, sim);
  return {
    params: fitted, values, n: obs.length,
    sse_before: before.sse, sse_after: after.sse, residuals: after.residuals
  };
}

/** Write a fitted parameter set as a loadable override module. */
export function writeOverrides(file, values) {
  const lines = Object.entries(values).map(
    ([k, v]) => '  ' + JSON.stringify(k) + ': ' + v.after.toPrecision(6)
      + ',   // was ' + Number(v.before).toPrecision(6) + '  (' + v.fold.toFixed(2) + 'x)'
      + (v.atBound ? '  <-- AT BOUND, not identified by the data' : '')
  );
  fs.writeFileSync(file,
    '// Generated by src/fit.js. Apply with applyOverrides(P, OVERRIDES).\n'
    + 'export const OVERRIDES = {\n' + lines.join('\n') + '\n};\n');
  return file;
}

/** Apply an override map (dotted paths) onto a parameter object in place. */
export function applyOverrides(par, overrides) {
  for (const [path, v] of Object.entries(overrides)) setPath(par, path, v);
  return par;
}
