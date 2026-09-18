import { rhs, observables } from '../src/model.js';
import { initialState, IDX } from '../src/state.js';
import { rk4 } from '../src/integrate.js';
import { P_GIWI, P_CHIR_PULSE, P_CHIR_CONTINUOUS, P_MAINTENANCE, protocol } from '../src/protocol.js';
import { writeTrajectory, table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// The central claim of the GiWi protocol (Lian et al. 2012 PNAS) is that
// cardiac differentiation needs canonical Wnt HIGH then LOW. Nothing in this
// model was fitted to that result: it should fall out of the single repression
// term rep(TCF_B) on NKX2-5, plus the GATA4 latch that carries cardiac-mesoderm
// identity across the gap. If it does not, the model is wrong.
// ---------------------------------------------------------------------------

const DAYS = 10, T = DAYS * 24;

const arms = {
  'no CHIR (control)':        P_MAINTENANCE,
  'CHIR 6uM d0-d1 only':      P_CHIR_PULSE(6, 24),
  'CHIR 6uM CONTINUOUS':      P_CHIR_CONTINUOUS(6),
  'GiWi: CHIR d0-1, IWP d3-5': P_GIWI({ dose: 6, chirEnd: 24, iwpStart: 72, iwpEnd: 120 }),
  'IWP d3-5 only (no CHIR)':  protocol([{ from: 72, to: 120, set: { iwp: 0.95 } }])
};

const results = {};
for (const [name, u] of Object.entries(arms)) {
  const traj = rk4(rhs, initialState(), { tEnd: T, h: 0.005, sample: 1, u });
  results[name] = traj;
}

// --- summary table ---------------------------------------------------------
const peak = (traj, key) => Math.max(...traj.map((p) => p.y[IDX[key]]));
const at   = (traj, key, h) => traj.find((p) => Math.abs(p.t - h) < 0.51).y[IDX[key]];

console.log('=== 03  Biphasic Wnt requirement for cardiac differentiation ===\n');
console.log(table(
  ['arm', 'peak TCF_B', 'peak TBXT', 'peak MESP1', 'GATA4 d10', 'NKX2-5 d10', 'TNNT2 d10', 'viab d10'],
  Object.entries(results).map(([n, tr]) => [
    n, r2(peak(tr, 'TCF_B')), r2(peak(tr, 'TBXT')), r2(peak(tr, 'MESP1')),
    r2(at(tr, 'GATA4', T)), r2(at(tr, 'NKX25', T)), r2(at(tr, 'TNNT2', T)), r2(at(tr, 'V', T))
  ])
));

// --- day-by-day trace for the GiWi arm ------------------------------------
console.log('\nGiWi arm, daily trace:');
const gi = results['GiWi: CHIR d0-1, IWP d3-5'];
console.log(table(
  ['day', 'GSK3act', 'Bcat_c', 'TCF_B', 'OCT4', 'SOX2', 'NANOG', 'TBXT', 'MESP1', 'GATA4', 'NKX2-5', 'TNNT2'],
  [...Array(DAYS + 1).keys()].map((dday) => {
    const p = gi.find((q) => Math.abs(q.t - dday * 24) < 0.51);
    const o = observables(p.y);
    return [dday, r2(o.GSK3_active), r2(p.y[IDX.Bcat_c], 1), r2(p.y[IDX.TCF_B]),
      r2(p.y[IDX.OCT4]), r2(p.y[IDX.SOX2]), r2(p.y[IDX.NANOG]), r2(p.y[IDX.TBXT]),
      r2(p.y[IDX.MESP1]), r2(p.y[IDX.GATA4]), r2(p.y[IDX.NKX25]), r2(p.y[IDX.TNNT2])];
  })
));

for (const [n, tr] of Object.entries(results)) {
  writeTrajectory(`03_${n.replace(/[^a-z0-9]+/gi, '_').toLowerCase()}.csv`, tr,
    (p) => observables(p.y));
}
console.log('\ntrajectories -> out/03_*.csv');
