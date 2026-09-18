import { P } from '../src/params.js';
import * as cardiac from '../src/model.js';
import { IDX, initialState } from '../src/state.js';
import { fit, makeSimulator, objective } from '../src/fit.js';
import { table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// Does the fitter actually recover parameters it is supposed to?
//
// Generate synthetic data from a KNOWN perturbed parameter set, hand only the
// data to the fitter, and see how close it gets back. If recovery is poor here
// it will be worse on real data, and the fix is more informative timepoints
// rather than more iterations.
// ---------------------------------------------------------------------------

const model = { rhs: cardiac.rhs, initialState, IDX };

const PROTOCOLS = {
  none:   [],
  chir2:  [{ from: 0, to: 24, set: { chir: 2 } }],
  chir4:  [{ from: 0, to: 24, set: { chir: 4 } }],
  chir6:  [{ from: 0, to: 24, set: { chir: 6 } }],
  chir12: [{ from: 0, to: 24, set: { chir: 12 } }],
  giwi6:  [{ from: 0, to: 24, set: { chir: 6 } }, { from: 72, to: 120, set: { iwp: 0.95 } }]
};

// --- the ground truth we will try to recover -------------------------------
const TRUTH = {
  'wnt.k_phos':     P.wnt.k_phos * 1.45,
  'wnt.k_syn_ax2':  P.wnt.k_syn_ax2 * 0.65,
  'chir.Kp':        P.chir.Kp * 1.30
};
const truthPar = JSON.parse(JSON.stringify(P));
for (const [path, v] of Object.entries(TRUTH)) {
  const k = path.split('.');
  truthPar[k[0]][k[1]] = v;
}

// --- synthetic observations: a TOPflash dose-response + a qPCR timecourse ---
const sim = makeSimulator(model, PROTOCOLS);
const obs = [];
for (const proto of ['chir2', 'chir4', 'chir6', 'chir12']) {
  obs.push({ kind: 'topflash', protocol: proto, readout: 'TCF_B', time_h: 24,
             value: sim.level(truthPar, proto, 'TCF_B', 24) / sim.level(truthPar, 'none', 'TCF_B', 24),
             sd: 0.20, ref: 'none' });
}
sim.reset();
for (const t of [6, 12, 24, 48]) {
  // Normalised to the undifferentiated control at the SAME timepoint, which is
  // what ddCt does. TBXT starts at zero, so a fold-change against t = 0 is
  // undefined - fit() will refuse that row rather than quietly returning Inf.
  obs.push({ kind: 'qpcr', protocol: 'giwi6', readout: 'TBXT', time_h: t,
             value: sim.level(truthPar, 'giwi6', 'TBXT', t) / sim.level(truthPar, 'none', 'TBXT', t),
             sd: 0.30, ref: 'none' });
  sim.reset();
}

console.log('=== 08  Fitter recovery check (synthetic data, known answer) ===\n');
console.log(`${obs.length} synthetic observations: 4 TOPflash doses + 4 qPCR timepoints\n`);

const t0 = Date.now();
const res = fit(model, PROTOCOLS, obs, Object.keys(TRUTH), { bounds: 4, restarts: 3 });
const secs = ((Date.now() - t0) / 1000).toFixed(1);

console.log(table(['parameter', 'start', 'truth', 'recovered', 'error', 'flag'],
  Object.entries(res.values).map(([k, v]) => {
    const err = 100 * (v.after - TRUTH[k]) / TRUTH[k];
    return [k, r2(v.before, 3), r2(TRUTH[k], 3), r2(v.after, 3),
      (err >= 0 ? '+' : '') + err.toFixed(1) + ' %',
      v.atBound ? 'AT BOUND' : (Math.abs(err) < 10 ? 'ok' : 'loose')];
  })));

console.log(`\nSSE  ${r2(res.sse_before, 3)}  ->  ${r2(res.sse_after, 4)}`
  + `   (${res.n} points, ${secs} s)`);

const maxErr = Math.max(...Object.entries(res.values)
  .map(([k, v]) => Math.abs(100 * (v.after - TRUTH[k]) / TRUTH[k])));
console.log(maxErr < 10
  ? `\nPASS - worst parameter recovered to ${maxErr.toFixed(1)} %. The fitter works;`
    + '\n       point it at data/ with your own numbers.'
  : `\nWEAK - worst parameter off by ${maxErr.toFixed(1)} %. These three are not`
    + '\n       separately identified by a dose-response plus one timecourse.'
    + '\n       Add timepoints in the first 12 h, where the three act differently.');
