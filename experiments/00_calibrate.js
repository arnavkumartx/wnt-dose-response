import { rhs, gsk3Activity, cellularIC50 } from '../src/model.js';
import { initialState, IDX, SPECIES } from '../src/state.js';
import { rk4, steadyState } from '../src/integrate.js';
import { P_MAINTENANCE, P_CHIR_PULSE, constant } from '../src/protocol.js';
import { table, r2 } from '../src/report.js';

console.log('=== 00  Calibration and numerical checks ===\n');

// --- 1. numerical convergence: halve the step, compare ---------------------
const u = P_CHIR_PULSE(6, 24);
const coarse = rk4(rhs, initialState(), { tEnd: 120, h: 0.01,   sample: 120, u });
const fine   = rk4(rhs, initialState(), { tEnd: 120, h: 0.0025, sample: 120, u });
const a = coarse[coarse.length - 1].y, b = fine[fine.length - 1].y;
let maxRel = 0, worst = '';
for (let i = 0; i < a.length; i++) {
  const scale = Math.max(Math.abs(b[i]), 1e-3);
  const rel = Math.abs(a[i] - b[i]) / scale;
  if (rel > maxRel) { maxRel = rel; worst = SPECIES[i]; }
}
console.log(`step halving h=0.01 -> 0.0025 : max relative drift ${(maxRel * 100).toExponential(2)} %  (${worst})`);
console.log(maxRel < 1e-4 ? '  PASS - h = 0.005 is converged\n' : '  WARN - tighten the step\n');

// --- 2. unperturbed hPSC steady state --------------------------------------
const ss = steadyState(rhs, initialState(), P_MAINTENANCE, { tEnd: 300 });
console.log('Undifferentiated steady state (no drug):');
console.log(table(['species', 'value', 'target'], [
  ['Bcat_c  (nM)', r2(ss[IDX.Bcat_c]), '~15-30'],
  ['Bcat_n  (nM)', r2(ss[IDX.Bcat_n]), 'low'],
  ['Bcat_m  (nM)', r2(ss[IDX.Bcat_m]), 'bulk of pool'],
  ['TCF_B   (nM)', r2(ss[IDX.TCF_B]),  '< 1.5'],
  ['AXIN2   (au)', r2(ss[IDX.AXIN2]),  'low'],
  ['OCT4', r2(ss[IDX.OCT4]), '> 0.85'],
  ['SOX2', r2(ss[IDX.SOX2]), '> 0.85'],
  ['NANOG', r2(ss[IDX.NANOG]), '> 0.70'],
  ['TBXT', r2(ss[IDX.TBXT]), '< 0.05'],
  ['TNNT2', r2(ss[IDX.TNNT2]), '~0']
]));

// --- 3. the Cheng-Prusoff shift --------------------------------------------
const { KiApp } = gsk3Activity(0);
console.log('\nGSK3 inhibition by CHIR99021:');
console.log(`  biochemical Ki (0 ATP)        ${(0.0099 * 1000).toFixed(1)} nM`);
console.log(`  ATP-corrected Ki_app (3 mM)   ${KiApp.toFixed(2)} uM   <- ${(KiApp / 0.0099).toFixed(0)}x shift`);
console.log(`  apparent MEDIUM IC50          ${cellularIC50().toFixed(2)} uM  (Kp = 0.55)`);
console.log(`  -> protocols use 3-12 uM because that is 1-6x the cellular IC50`);

// --- 4. fold-induction of the Wnt reporter ---------------------------------
console.log('\nSteady-state TCF:beta-catenin vs sustained CHIR (TOPflash analogue):');
const base = steadyState(rhs, initialState(), constant({ chir: 0 }), { tEnd: 300 })[IDX.TCF_B];
const rows = [];
for (const dose of [0, 0.5, 1, 2, 3, 6, 9, 12]) {
  const s = steadyState(rhs, initialState(), constant({ chir: dose }), { tEnd: 300 });
  const g = gsk3Activity(s[IDX.CHIR_in]).active;
  rows.push([`${dose} uM`, r2(100 * (1 - g), 1) + ' %', r2(s[IDX.TCF_B]),
             r2(s[IDX.TCF_B] / base, 1) + 'x', r2(s[IDX.Bcat_c]), r2(s[IDX.V], 2)]);
}
console.log(table(['CHIR', 'GSK3 inhib', 'TCF_B nM', 'fold', 'Bcat_c nM', 'viab'], rows));
