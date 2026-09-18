import * as H from '../src/hsc.js';
import { S } from '../src/wnt-core.js';
import { rk4 } from '../src/integrate.js';
import { protocol } from '../src/protocol.js';
import { writeCSV, table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// Cord blood CD34+, 7-day ex vivo culture with SCF/TPO/FLT3L +/- CHIR99021.
//
// What the literature says should happen (PubMed):
//   Holmes 2008  10.1634/stemcells.2007-0600  GSK-3b inhibition DELAYS CD34+
//                expansion but PRESERVES repopulating activity; CXCR4 up.
//   Ko 2011      10.1002/stem.551             lengthens cycle time, p57 up,
//                cyclin D1 down, better engraftment per repopulating unit.
//   Luis 2011    10.1016/j.stem.2011.07.017   the Wnt dose-response is
//                non-monotonic; push too hard and you deplete the HSC pool.
//
// So the endpoint that matters is not "how much did it expand" but
// "how many engrafting units came out per input cell".
// ---------------------------------------------------------------------------
const T = 7 * 24;
const DOSES = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12];

const run = (dose) => {
  const u = protocol([{ from: 0, to: T, set: { chir: dose, cytokines: 1.0 } }]);
  const tr = rk4(H.rhs, H.initialState(), { tEnd: T, h: 0.01, sample: 24, u });
  return { tr, y: tr[tr.length - 1].y };
};

console.log('=== 06  Cord blood CD34+, 7-day culture, CHIR dose-response ===\n');

const rows = [], csv = [];
for (const dose of DOSES) {
  const { y } = run(dose);
  const e = H.engraftment(y);
  rows.push([`${dose}`, r2(100 * (1 - e.GSK3_active), 0), r2(y[S.TCF_B]),
    r2(y[H.IDX.HOXB4]), r2(y[H.IDX.MYC]), r2(y[H.IDX.CDKN1C]), r2(y[H.IDX.CCND1]),
    r2(e.fold_expansion, 1), r2(e.stem_frequency, 3), r2(e.homing), r2(e.engraftment, 2)]);
  csv.push([dose, y[S.TCF_B], y[H.IDX.HOXB4], y[H.IDX.MYC], y[H.IDX.CDKN1C], y[H.IDX.CCND1],
    y[H.IDX.CXCR4], y[H.IDX.CD34], y[H.IDX.CD90], y[H.IDX.CD45RA], y[H.IDX.LYMPH],
    e.fold_expansion, e.stem_frequency, e.homing, e.units, e.engraftment]);
}
console.log(table(['CHIR uM', '%GSK3i', 'TCF_B', 'HOXB4', 'MYC', 'p57', 'CCND1',
  'fold exp', 'STEM', 'homing', 'ENGRAFT'], rows));

const best = csv.reduce((a, b) => (b[15] > a[15] ? b : a));
const ctrl = csv[0];
console.log(`\ncontrol (0 uM)   fold ${r2(ctrl[11], 1)}x,  STEM ${r2(ctrl[12], 3)},  engraftment ${r2(ctrl[15], 2)}`);
console.log(`optimum (${best[0]} uM)  fold ${r2(best[11], 1)}x,  STEM ${r2(best[12], 3)},  engraftment ${r2(best[15], 2)}`
  + `   -> ${r2(best[15] / ctrl[15], 2)}x control`);
console.log(`high dose (12 uM) fold ${r2(csv[12][11], 1)}x,  STEM ${r2(csv[12][12], 3)},  engraftment ${r2(csv[12][15], 2)}`
  + `   -> ${r2(csv[12][15] / ctrl[15], 2)}x control`);

console.log('\n  The optimum is INTERIOR and it is low. More CHIR buys more');
console.log('  self-renewal signalling and fewer cells, and past the MYC threshold');
console.log('  it buys neither. High-dose CHIR is worse than no CHIR at all.');

// --- lineage output --------------------------------------------------------
console.log('\nLineage output after 7 days:');
console.log(table(['CHIR uM', 'CD34', 'CD90', 'CD45RA', 'GATA1 (ery)', 'SPI1 (mye)', 'LYMPH'],
  csv.filter((r) => [0, 1.5, 3, 6, 12].includes(r[0])).map((r) => {
    const { y } = run(r[0]);
    return [r2(r[0], 1), r2(r[7]), r2(r[8]), r2(r[9]),
      r2(y[H.IDX.GATA1]), r2(y[H.IDX.SPI1]), r2(r[10])];
  })));
console.log('\n  Lymphoid output falls monotonically with Wnt dose - the model has');
console.log('  no lymphoid-friendly CHIR window at all, matching the persistent');
console.log('  T-cell depletion reported in Shen 2014.');

writeCSV('06_cordblood_dose.csv', csv,
  ['chir_uM', 'TCF_B', 'HOXB4', 'MYC', 'CDKN1C', 'CCND1', 'CXCR4', 'CD34', 'CD90',
   'CD45RA', 'LYMPH', 'fold_expansion', 'stem_frequency', 'homing', 'units', 'engraftment']);
console.log('\n-> out/06_cordblood_dose.csv');
