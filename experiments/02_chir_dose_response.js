import { runPopulation } from '../src/population.js';
import { P_GIWI, P_CHIR_PULSE } from '../src/protocol.js';
import { writeCSV, table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// Dose x duration response surface, scored the way a flow cytometer scores it:
// % of surviving cells above a marker threshold, over a heterogeneous dish.
// ---------------------------------------------------------------------------
const N = 150, T = 240, H = 0.02;
const run = (u, seed) => runPopulation(u, { n: N, tEnd: T, h: H, seed }).summary;

console.log('=== 02  CHIR dose-response (population of %d cells/condition) ==='.replace('%d', N) + '\n');

// --- 1-D sweep at 24 h -----------------------------------------------------
const DOSES = [0, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25];
const bar = (x) => '#'.repeat(Math.round(x / 3));
const sweep = DOSES.map((dose) => {
  const s = run(P_GIWI({ dose, chirEnd: 24 }), 1234);
  return { dose, ...s };
});

console.log('GiWi schedule, CHIR 0-24 h, IWP2 72-120 h, scored at day 10:\n');
console.log(table(['CHIR uM', '%TBXT', '%GATA4', '%NKX2-5', '%TNNT2+', 'viab', 'yield/input', ''],
  sweep.map((s) => [r2(s.dose, 1), r2(s.pct_TBXT, 0), r2(s.pct_GATA4, 0), r2(s.pct_NKX25, 0),
    r2(s.pct_TNNT2, 1), r2(s.viability, 2), r2(s.yield_TNNT2, 1), bar(s.yield_TNNT2)])));

const best = sweep.reduce((a, b) => (b.yield_TNNT2 > a.yield_TNNT2 ? b : a));
const window = sweep.filter((s) => s.yield_TNNT2 > 0.5 * best.yield_TNNT2).map((s) => s.dose);
console.log(`\noptimum ${best.dose} uM (${r2(best.yield_TNNT2, 1)} % yield);`
  + ` >50%-of-max window = ${Math.min(...window)}-${Math.max(...window)} uM`);

// --- dose x duration grid --------------------------------------------------
console.log('\n% TNNT2+ over dose x CHIR duration:\n');
const DURS = [12, 18, 24, 36, 48, 72];
const gridDoses = [1.5, 2, 3, 4, 6, 8, 12, 16, 20];
const grid = [];
const rows = gridDoses.map((dose) => [`${dose} uM`, ...DURS.map((durH) => {
  const s = run(P_GIWI({ dose, chirEnd: durH, iwpStart: Math.max(durH + 24, 72),
                         iwpEnd: Math.max(durH + 24, 72) + 48 }), 1234);
  grid.push([dose, durH, s.pct_TBXT, s.pct_GATA4, s.pct_NKX25, s.pct_TNNT2, s.viability, s.yield_TNNT2]);
  return r2(s.pct_TNNT2, 0);
})]);
console.log(table(['CHIR', ...DURS.map((h) => `${h}h`)], rows));
console.log('\n  The ridge runs diagonally: a shorter pulse needs a higher dose.');
console.log('  What the switch integrates is exposure above threshold, not dose alone.');

// --- how much of the window is set by heterogeneity alone? -----------------
console.log('\nEffect of clonal heterogeneity on the width of the usable window:');
const cvRows = [];
for (const scale of [0.25, 0.5, 1, 2]) {
  const cv = Object.fromEntries(Object.entries(
    { 'chir.Kp': 0.25, 'wnt.k_phos': 0.30, 'wnt.k_syn_b': 0.15, 'wnt.k_sec_wnt': 0.35, 'grn.k_tbxt': 0.20 }
  ).map(([k, v]) => [k, v * scale]));
  const pts = [1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 20].map((dose) => ({
    dose, y: runPopulation(P_GIWI({ dose, chirEnd: 24 }), { n: N, tEnd: T, h: H, seed: 99, cv }).summary.yield_TNNT2
  }));
  const mx = Math.max(...pts.map((p) => p.y));
  const w = pts.filter((p) => p.y > 0.5 * mx).map((p) => p.dose);
  cvRows.push([`${(scale * 100).toFixed(0)} % of default CV`, r2(mx, 1),
    `${Math.min(...w)}-${Math.max(...w)} uM`, r2(Math.max(...w) / Math.min(...w), 1) + 'x']);
}
console.log(table(['heterogeneity', 'peak yield %', 'usable window', 'fold-width'], cvRows));
console.log('\n  PREDICTION: a tighter clone has a HIGHER peak but a NARROWER window.');
console.log('  Lines that "work at any dose" are the heterogeneous ones, and they');
console.log('  cap out lower. That trade-off is testable by single-cell readout.');

writeCSV('02_dose_response_24h.csv',
  sweep.map((s) => [s.dose, s.pct_TBXT, s.pct_GATA4, s.pct_NKX25, s.pct_TNNT2, s.viability, s.yield_TNNT2]),
  ['chir_uM', 'pct_TBXT', 'pct_GATA4', 'pct_NKX25', 'pct_TNNT2', 'viability', 'yield_per_input_pct']);
writeCSV('02_dose_duration_grid.csv', grid,
  ['chir_uM', 'duration_h', 'pct_TBXT', 'pct_GATA4', 'pct_NKX25', 'pct_TNNT2', 'viability', 'yield_pct']);
console.log('\n-> out/02_dose_response_24h.csv, out/02_dose_duration_grid.csv');
