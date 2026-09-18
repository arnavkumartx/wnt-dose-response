import { runPopulation } from '../src/population.js';
import { protocol } from '../src/protocol.js';
import { writeCSV, table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// Wnt alone does not pick a lineage. In the model, TBXT reads Wnt, but EOMES
// is an AND-gate on Wnt *and* SMAD2/3, and SOX17 sits downstream of EOMES.
// So the mesoderm/endoderm split should be governed by the Activin axis at
// constant CHIR - which is exactly how the two protocols differ in practice.
// ---------------------------------------------------------------------------
// TIMING MATTERS MORE THAN IT LOOKS. MESP1 is a transient pulse: it peaks at
// 24-26 h, is down to ~0.3 by 48 h, and is gone by day 3. Scored at day 7 it
// reads as the leak floor at EVERY dose and says nothing at all; scored at 48 h
// it is already flat. TBXT peaks at 24 h, so 26 h catches both. Each marker is
// read when it is actually informative, which for a transient gene is doing as
// much work as the dose axis.
const N = 120, H = 0.02;
const T_STREAK = 26, T_LATE = 168;

const CHIRS   = [0, 2, 3, 4, 6, 9];
const ACTIVIN = [0.18, 0.3, 0.5, 0.75, 1.0];   // 0.18 = autocrine floor, 1.0 = 100 ng/mL

console.log('=== 05  Wnt x Nodal lineage map ===\n');
console.log('CHIR 0-24 h; Activin A held for 72 h.');
console.log('TBXT and MESP1 scored at 26 h, where both peak; SOX17 at day 7.\n');

const grid = [];
// NOTE ON READOUTS: GATA4 is deliberately NOT used to score the mesoderm side.
// GATA4/GATA6 are expressed in definitive endoderm as well as cardiac mesoderm,
// so GATA4+ does not discriminate the two lineages - in this model it goes high
// in both corners. MESP1 is the cardiac-mesoderm-specific call.
const markers = ['TBXT', 'SOX17', 'MESP1', 'GATA4', 'NKX25', 'OCT4'];
const makeProto = (dose, act) => protocol([
  { from: 0, to: 24, set: { chir: dose, activin: act } },
  { from: 24, to: 72, set: { activin: act } }
]);

// One population run per (dose, Activin, timepoint), cached so each panel reads
// the same runs instead of repeating them.
const cache = new Map();
const summary = (dose, act, tEnd) => {
  const key = `${dose}|${act}|${tEnd}`;
  if (!cache.has(key)) {
    cache.set(key, runPopulation(makeProto(dose, act),
      { n: N, tEnd, h: H, seed: 77, markers }).summary);
  }
  return cache.get(key);
};

const PANELS = [
  ['TBXT  (streak, 26 h)',               'TBXT',  T_STREAK],
  ['MESP1 (cardiac mesoderm, 26 h)',     'MESP1', T_STREAK],
  ['SOX17 (definitive endoderm, day 7)', 'SOX17', T_LATE]
];

// Both % positive and mean level, because MESP1 peaks at ~0.57 against a 0.5
// gate: the percentage sits on a knife edge there and swings on small parameter
// changes, while the mean does not. Reporting only the percentage would make a
// gating artefact look like a biological effect.
for (const [name, marker, tEnd] of PANELS) {
  for (const [stat, prefix, dp] of [['% positive', 'pct_', 0], ['mean level', 'mean_', 2]]) {
    console.log(`${name} - ${stat} :`);
    console.log(table(['', ...ACTIVIN.map((a) => `Act ${a}`)],
      CHIRS.map((dose) => [`CHIR ${dose}`,
        ...ACTIVIN.map((act) => r2(summary(dose, act, tEnd)[prefix + marker], dp))])));
    console.log('');
  }
}

for (const dose of CHIRS) for (const act of ACTIVIN) {
  const e = summary(dose, act, T_STREAK), l = summary(dose, act, T_LATE);
  grid.push([dose, act, e.pct_TBXT, e.pct_MESP1, l.pct_SOX17, l.pct_GATA4,
             l.pct_NKX25, l.viability]);
}

console.log('  Reading: Wnt (rows) sets whether cells leave pluripotency at all.');
console.log('  Nodal (columns) sets which side of the streak they leave through.');
console.log('  Neither knob alone picks a lineage - EOMES is an AND-gate on both.');
console.log('\n  Both arms need a LATCH to survive the CHIR washout: GATA4 on the');
console.log('  cardiac side, SOX17/FOXA2 on the endoderm side. Remove either and');
console.log('  that lineage collapses when the 24 h spike ends, on any schedule.');

writeCSV('05_lineage_map.csv', grid,
  ['chir_uM', 'activin', 'pct_TBXT_26h', 'pct_MESP1_26h', 'pct_SOX17_d7',
   'pct_GATA4_d7', 'pct_NKX25_d7', 'viability_d7']);
console.log('\n-> out/05_lineage_map.csv');
