import { P } from '../src/params.js';
import { rhs } from '../src/model.js';
import { initialState, IDX } from '../src/state.js';
import { rk4 } from '../src/integrate.js';
import { P_GIWI } from '../src/protocol.js';
import { writeCSV, table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// Local sensitivity: perturb each parameter +/-20 % and measure the effect on
// two endpoints, to identify which [A]-tagged assumptions matter enough to be
// worth measuring.
// ---------------------------------------------------------------------------
const T = 240;
const clone = (o) => JSON.parse(JSON.stringify(o));

function endpoints(par) {
  const tr = rk4(rhs, initialState(), { tEnd: T, h: 0.02, sample: 24, u: P_GIWI({ dose: 6 }), par });
  const last = tr[tr.length - 1].y;
  const peakTCF = Math.max(...tr.map((p) => p.y[IDX.TCF_B]));
  return { TNNT2: last[IDX.TNNT2], peakTCF };
}

const base = endpoints(P);
console.log('=== 04  Local sensitivity (+/-20 %), GiWi 6 uM ===\n');
console.log(`baseline: TNNT2(d10) = ${r2(base.TNNT2)},  peak TCF_B = ${r2(base.peakTCF)} nM\n`);

const targets = [];
for (const grp of ['chir', 'wnt', 'grn']) {
  for (const k of Object.keys(P[grp])) {
    if (typeof P[grp][k] === 'number' && P[grp][k] !== 0) targets.push([grp, k]);
  }
}

const rows = [];
for (const [grp, k] of targets) {
  const s = {};
  for (const f of [0.8, 1.2]) {
    const par = clone(P);
    par[grp][k] = P[grp][k] * f;
    // GRN synthesis/degradation pairs must stay matched to keep TFs on [0,1]
    if (grp === 'grn') {
      const twin = (k.startsWith('k_') ? 'd_' : k.startsWith('d_') ? 'k_' : null);
      if (twin && (twin + k.slice(2)) in par.grn) par.grn[twin + k.slice(2)] = par[grp][k];
    }
    s[f] = endpoints(par);
  }
  // normalised sensitivity: (dY/Y) / (dX/X) over the +/-20 % span
  const sTNNT2 = ((s[1.2].TNNT2 - s[0.8].TNNT2) / Math.max(base.TNNT2, 1e-6)) / 0.4;
  const sTCF   = ((s[1.2].peakTCF - s[0.8].peakTCF) / Math.max(base.peakTCF, 1e-6)) / 0.4;
  rows.push({ name: `${grp}.${k}`, sTNNT2, sTCF, mag: Math.abs(sTNNT2) + Math.abs(sTCF) });
}

rows.sort((a, b) => b.mag - a.mag);
console.log('Top 18 parameters by normalised sensitivity  ( (dY/Y)/(dX/X) ):\n');
console.log(table(['parameter', 'S[TNNT2 d10]', 'S[peak TCF_B]'],
  rows.slice(0, 18).map((r) => [r.name, r2(r.sTNNT2), r2(r.sTCF)])));

const stiff = rows.filter((r) => Math.abs(r.sTNNT2) < 0.02 && Math.abs(r.sTCF) < 0.02);
console.log(`\n${stiff.length}/${rows.length} parameters move neither endpoint by >2 % per 10 % change.`);
console.log('Those are safe to leave at their assumed values; the list above is');
console.log('where measurement effort should go.');

writeCSV('04_sensitivity.csv', rows.map((r) => [r.name, r.sTNNT2, r.sTCF]),
  ['parameter', 'S_TNNT2_d10', 'S_peak_TCF_B']);
console.log('\n-> out/04_sensitivity.csv');
