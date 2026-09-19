import { gsk3Activity, cellularIC50 } from '../src/model.js';
import { P } from '../src/params.js';
import { writeCSV, table, r2 } from '../src/report.js';

// ---------------------------------------------------------------------------
// The gap between the biochemical IC50 of CHIR99021 and the dose used at the
// bench.
// ---------------------------------------------------------------------------
console.log('=== 01  CHIR99021 -> GSK3 occupancy ===\n');

const { KiApp } = gsk3Activity(0);
console.log(`Ki (biochemical, ATP-free)   ${(P.chir.Ki_gsk3b * 1000).toFixed(1)} nM`);
console.log(`Km(ATP) for GSK3-beta        ${P.chir.Km_ATP} uM`);
console.log(`cellular [ATP]               ${P.chir.ATP_cell / 1000} mM`);
console.log(`Ki_app = Ki*(1+ATP/Km)       ${KiApp.toFixed(2)} uM    (${(KiApp / P.chir.Ki_gsk3b).toFixed(0)}x shift)`);
console.log(`cell:medium partition Kp     ${P.chir.Kp}`);
console.log(`=> MEDIUM IC50               ${cellularIC50().toFixed(2)} uM\n`);

const rows = [];
for (const dose of [0, 0.3, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20]) {
  const cIn = P.chir.Kp * dose;
  const g = gsk3Activity(cIn).active;
  rows.push([`${dose}`, r2(cIn), r2(100 * (1 - g), 1), r2(g, 3)]);
}
console.log(table(['CHIR medium uM', 'intracell uM', '% GSK3 inhibited', 'GSK3 activity'], rows));

// ATP sensitivity: why the same dose behaves differently between batches
console.log('\nSensitivity of the effective dose to the cellular ATP pool:');
const atpRows = [];
const saved = P.chir.ATP_cell;
for (const atp of [1000, 2000, 3000, 4000, 5000]) {
  P.chir.ATP_cell = atp;
  const need = gsk3Activity(0).KiApp / P.chir.Kp;   // dose for 50% inhibition
  const at6 = gsk3Activity(P.chir.Kp * 6).active;
  atpRows.push([`${atp / 1000} mM`, r2(need), r2(100 * (1 - at6), 1) + ' %']);
}
P.chir.ATP_cell = saved;
console.log(table(['cellular ATP', 'medium IC50 uM', 'inhibition at 6 uM'], atpRows));
console.log('\n  Metabolic state shifts the effective CHIR dose ~3-fold across a');
console.log('  physiological ATP range. Testable: clamp ATP (2-DG, oligomycin,');
console.log('  galactose medium) and the CHIR dose-response should slide.');

const sweep = [];
for (let d = 0; d <= 20; d += 0.1) {
  const g = gsk3Activity(P.chir.Kp * d).active;
  sweep.push([d, P.chir.Kp * d, g, 100 * (1 - g)]);
}
writeCSV('01_gsk3_occupancy.csv', sweep, ['chir_medium_uM', 'chir_intracell_uM', 'gsk3_activity', 'pct_inhibited']);
console.log('\nsweep -> out/01_gsk3_occupancy.csv');
