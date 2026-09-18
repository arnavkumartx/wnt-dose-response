import { rhs } from '../src/model.js';
import { initialState, IDX } from '../src/state.js';
import { S } from '../src/wnt-core.js';
import { rk4 } from '../src/integrate.js';
import { protocol } from '../src/protocol.js';
import { destructionComplex, nuclearOutput } from '../src/interfaces.js';
import { table, r2, writeCSV } from '../src/report.js';

// ---------------------------------------------------------------------------
// A prediction that only exists once the destruction complex is explicit.
//
// IWP2  (porcupine)  blocks WNT palmitoylation -> no new ligand. It acts ABOVE
//                    the receptor, so it cannot touch CHIR at all.
// IWR-1 (tankyrase)  blocks AXIN1 PARsylation -> Axin is degraded less, total
//                    scaffold rises. More scaffold means more residual complex
//                    even when each bound GSK3 is mostly inhibited, so it acts
//                    BELOW CHIR and partially reverses it.
//
// The two are used almost interchangeably in cardiac protocols. The model says
// they are not interchangeable if CHIR is still present.
// ---------------------------------------------------------------------------
console.log('=== 07  IWR-1 and IWP2 are not interchangeable ===\n');

console.log('Destruction-complex activity (a.u.), scaffold 1.2:');
console.log(table(['condition', 'Axin', 'assembled', 'ACTIVE', 'vs untreated'],
  [['untreated', 1.0, 0], ['CHIR 6 uM', 0.376, 0], ['CHIR 6 uM + IWR-1', 0.376, 0.90],
   ['IWR-1 alone', 1.0, 0.90]].map(([lab, g, iwr]) => {
    const d = destructionComplex(1.2, g, iwr);
    const base = destructionComplex(1.2, 1.0, 0).active;
    return [lab, r2(d.axin), r2(d.assembled, 3), r2(d.active, 3), r2(d.active / base, 2) + 'x'];
  })));

// --- does it change the fate outcome? --------------------------------------
const T = 10 * 24;
const arms = {
  'CHIR d0-1, nothing after':        [],
  'CHIR d0-1, IWP2 d3-5':            [{ from: 72, to: 120, set: { iwp: 0.95 } }],
  'CHIR d0-1, IWR-1 d3-5':           [{ from: 72, to: 120, set: { iwr: 0.90 } }],
  'CHIR CONTINUOUS + IWP2 d3-5':     [{ from: 0, to: T, set: { chir: 6 } }, { from: 72, to: 120, set: { chir: 6, iwp: 0.95 } }],
  'CHIR CONTINUOUS + IWR-1 d3-5':    [{ from: 0, to: T, set: { chir: 6 } }, { from: 72, to: 120, set: { chir: 6, iwr: 0.90 } }]
};

const rows = [], csv = [];
for (const [name, extra] of Object.entries(arms)) {
  const base = name.startsWith('CHIR CONTINUOUS') ? [] : [{ from: 0, to: 24, set: { chir: 6 } }];
  const u = protocol([...base, ...extra]);
  const tr = rk4(rhs, initialState(), { tEnd: T, h: 0.005, sample: 6, u });
  const y = tr[tr.length - 1].y;
  const minW = Math.min(...tr.filter((p) => p.t >= 72 && p.t <= 120).map((p) => p.y[S.TCF_B]));
  rows.push([name, r2(minW), r2(y[IDX.GATA4]), r2(y[IDX.NKX25]), r2(y[IDX.TNNT2])]);
  csv.push([name, minW, y[IDX.GATA4], y[IDX.NKX25], y[IDX.TNNT2]]);
}
console.log('\nCardiac outcome at day 10:');
console.log(table(['arm', 'min TCF_B in window', 'GATA4', 'NKX2-5', 'TNNT2'], rows));

console.log('\n  With CHIR withdrawn, the two inhibitors do the same job.');
console.log('  With CHIR still on at 6 uM, IWP2 does nothing at all (it acts above');
console.log('  the receptor); IWR-1 lowers TCF_B by ~29 % but NOT below the NKX2-5');
console.log('  threshold, so neither rescues. The question is where the crossover is.');

// --- where CAN a tankyrase inhibitor rescue continuous CHIR? ---------------
console.log('\nContinuous CHIR at varying dose, inhibitor d3-5. TNNT2 at day 10:');
const sweepRows = [1, 1.5, 2, 3, 4, 6, 9].map((dose) => {
  const row = [String(dose)];
  for (const extra of [{}, { iwp: 0.95 }, { iwr: 0.90 }]) {
    const u = protocol([
      { from: 0, to: T, set: { chir: dose } },
      { from: 72, to: 120, set: { chir: dose, ...extra } }
    ]);
    const tr = rk4(rhs, initialState(), { tEnd: T, h: 0.005, sample: T, u });
    row.push(r2(tr[tr.length - 1].y[IDX.TNNT2]));
  }
  return row;
});
console.log(table(['CHIR uM', 'no inhibitor', 'IWP2', 'IWR-1'], sweepRows));
console.log('\n  There is NO rescue window. A 2.35x scaffold increase cannot outrun');
console.log('  direct GSK3 inhibition, and below ~2 uM the streak never fires in');
console.log('  the first place. So the explicit-scaffold layer does NOT predict');
console.log('  the two inhibitors give different fates in the standard protocol.');
console.log('\n  Where it does separate them is the DEPTH of the Wnt drop during');
console.log('  the washout window: IWR-1 reaches TCF_B 0.82 against IWP2 1.55, a');
console.log('  ~2x difference, because it suppresses below the basal set point');
console.log('  rather than only removing ligand.');
console.log('\n  TESTABLE, and cheap: run the normal GiWi protocol with each');
console.log('  inhibitor and read AXIN2 by qPCR (or TOPflash) at d4. The model');
console.log('  says IWR-1 drives it roughly twice as far down as IWP2 while TNNT2');
console.log('  comes out the same. If the two are indistinguishable on AXIN2,');
console.log('  Axin is not limiting in these cells and the scaffold term is wrong.');

// --- the other interface: ICAT titration ------------------------------------
console.log('\nARM-groove competition at 12 nM nuclear beta-catenin:');
console.log(table(['ICAT nM', 'free bcat', 'bcat:TCF', 'W_eff', '% of no-ICAT'],
  [0, 2, 5, 10, 20, 40].map((icat) => {
    const r = nuclearOutput(12, { ICAT_tot: icat });
    const b = nuclearOutput(12, { ICAT_tot: 0 });
    return [icat, r2(r.Bfree), r2(r.bcatTCF), r2(r.W_eff), r2(100 * r.W_eff / b.W_eff, 0)];
  })));
console.log('\n  ICAT and TCF are mutually exclusive on ARM 3-9, so ICAT titrates');
console.log('  the output down without changing beta-catenin levels at all - a');
console.log('  western blot would look unchanged while the reporter falls.');

writeCSV('07_iwr_vs_iwp.csv', csv, ['arm', 'min_TCF_B', 'GATA4', 'NKX25', 'TNNT2']);
console.log('\n-> out/07_iwr_vs_iwp.csv');
