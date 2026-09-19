#!/usr/bin/env node
// Run a protocol without editing any source.
//
//   node run.js --chir 6 --chir-end 24 --iwp 72:120 --days 10
//   node run.js --chir 4 --chir-end 48 --activin 1.0 --days 7      # endoderm arm
//   node run.js --chir 6 --chir-end 24 --cells 200                 # dish, not one cell
//   node run.js --model hsc --chir 1 --chir-end 168 --days 7       # cord blood CD34+
//   node run.js --chir 6 --chir-end 24 --csv mine.csv

import * as cardiac from './src/model.js';
import * as cardiacState from './src/state.js';
import * as hspc from './src/hsc.js';
import { gsk3Activity } from './src/wnt-core.js';
import { rk4 } from './src/integrate.js';
import { protocol } from './src/protocol.js';
import { runPopulation } from './src/population.js';
import { writeTrajectory, table, r2 } from './src/report.js';

const argv = process.argv.slice(2);
const arg = (name, dflt) => {
  const i = argv.indexOf(`--${name}`);
  return i >= 0 && argv[i + 1] && !argv[i + 1].startsWith('--') ? argv[i + 1] : dflt;
};
const num = (name, dflt) => Number(arg(name, dflt));
const has = (name) => argv.includes(`--${name}`);

if (has('help') || has('h')) {
  console.log(`
  --chir <uM>          CHIR99021 dose                       (default 6)
  --chir-start <h>     when CHIR goes on                    (default 0)
  --chir-end <h>       when CHIR comes off                  (default 24)
  --iwp <from:to>      Wnt-inhibitor window in hours        (default 72:120, "none" to omit)
  --wnt3a <nM>         recombinant WNT3A instead of CHIR    (default 0)
  --activin <0..1>     Activin A / Nodal (1.0 = 100 ng/mL)  (default 0.18 = autocrine floor)
  --bmp4 <0..1>        BMP4                                 (default 0.10)
  --days <n>           length of the run                    (default 10)
  --cells <n>          ensemble size; 1 = deterministic      (default 1)
  --model <name>       cardiac (default) or hsc (cord blood CD34+)
  --cytokines <0..1>   SCF/TPO/FLT3L, hsc model only         (default 1.0)
  --csv <file>         write the trajectory to out/<file>
`);
  process.exit(0);
}

const modelName = arg('model', 'cardiac');
if (!['cardiac', 'hsc'].includes(modelName)) {
  console.error(`unknown --model "${modelName}"; expected cardiac or hsc`);
  process.exit(1);
}
const isHSC = modelName === 'hsc';
const M = isHSC
  ? { rhs: hspc.rhs, initialState: hspc.initialState, IDX: hspc.IDX, observables: null }
  : { rhs: cardiac.rhs, initialState: cardiacState.initialState, IDX: cardiacState.IDX,
      observables: cardiac.observables };

const dose      = num('chir', 6);
const chirStart = num('chir-start', 0);
const chirEnd   = num('chir-end', 24);
const iwpArg    = arg('iwp', '72:120');
const activin   = num('activin', 0.18);
const bmp4      = num('bmp4', 0.10);
const wnt3a     = num('wnt3a', 0);
const days      = num('days', 10);
const cells     = Math.max(1, Math.round(num('cells', 1)));
const T = days * 24;

const cytokines = num('cytokines', 1.0);

const windows = [];
const base = isHSC ? { cytokines } : { activin, bmp4 };
if (dose > 0 || wnt3a > 0) windows.push({ from: chirStart, to: chirEnd, set: { chir: dose, wnt3a, ...base } });
if (!isHSC && (activin !== 0.18 || bmp4 !== 0.10)) windows.push({ from: chirEnd, to: T, set: base });
if (!isHSC && iwpArg !== 'none') {
  const [a, b] = iwpArg.split(':').map(Number);
  windows.push({ from: a, to: b, set: { iwp: 0.95, ...base } });
}
const u = protocol(windows);

const SHOW = isHSC
  ? ['Bcat_c', 'TCF_B', 'HOXB4', 'MYC', 'CDKN1C', 'CCND1', 'CXCR4',
     'CD34', 'CD90', 'CD45RA', 'LYMPH', 'STEM', 'Ncell', 'V']
  : ['Bcat_c', 'TCF_B', 'WNT_s', 'OCT4', 'SOX2', 'NANOG', 'TBXT',
     'MESP1', 'SOX17', 'GATA4', 'NKX25', 'TNNT2', 'V'];
const MARKERS = isHSC
  ? ['CD34', 'CD90', 'CD45RA', 'LYMPH']
  : ['TBXT', 'GATA4', 'NKX25', 'TNNT2', 'SOX17', 'OCT4'];

// --json emits a machine-readable blob on stdout instead of tables, so the
// Python wrapper in notebooks/ can parse it directly rather than scraping
// formatted output. Nothing else prints in this mode.
const asJson = has('json');
const SPECIES_ALL = isHSC ? hspc.SPECIES : cardiacState.SPECIES;
const protocolSummary = {
  model: modelName, chir: dose, chir_start: chirStart, chir_end: chirEnd,
  wnt3a, days, cells, windows,
  ...(isHSC ? { cytokines } : { activin, bmp4, iwp: iwpArg })
};

if (asJson) {
  if (cells === 1) {
    const tr = rk4(M.rhs, M.initialState(), { tEnd: T, h: 0.005, sample: 1, u });
    const out = { mode: 'trajectory', protocol: protocolSummary, t_h: [], GSK3_active: {} };
    const cols = {};
    for (const s of SPECIES_ALL) cols[s] = [];
    const gsk = [];
    for (const p of tr) {
      out.t_h.push(p.t);
      gsk.push(gsk3Activity(p.y[M.IDX.CHIR_in]).active);
      for (const s of SPECIES_ALL) cols[s].push(p.y[M.IDX[s]]);
    }
    out.GSK3_active = gsk;
    out.species = cols;
    if (isHSC) out.engraftment = hspc.engraftment(tr[tr.length - 1].y);
    process.stdout.write(JSON.stringify(out));
  } else {
    const { summary } = runPopulation(u, { n: cells, tEnd: T, h: 0.02, model: M, markers: MARKERS });
    process.stdout.write(JSON.stringify({ mode: 'population', protocol: protocolSummary, summary }));
  }
  process.exit(0);
}

console.log(`\nmodel: ${modelName}`);
console.log(`protocol: CHIR ${dose} uM ${chirStart}-${chirEnd} h`
  + (wnt3a ? `, WNT3A ${wnt3a} nM` : '')
  + (isHSC ? `, cytokines ${cytokines}`
           : `, Activin ${activin}, BMP4 ${bmp4}`
             + (iwpArg === 'none' ? ', no Wnt inhibitor' : `, IWP ${iwpArg} h`))
  + `, ${days} days, ${cells} cell${cells > 1 ? 's' : ''}\n`);

if (cells === 1) {
  const tr = rk4(M.rhs, M.initialState(), { tEnd: T, h: 0.005, sample: 1, u });
  console.log(table(['day', 'GSK3act', ...SHOW],
    [...Array(days + 1).keys()].map((d) => {
      const p = tr.find((q) => Math.abs(q.t - d * 24) < 0.51);
      const act = gsk3Activity(p.y[M.IDX.CHIR_in]).active;
      return [d, r2(act), ...SHOW.map((s) => r2(p.y[M.IDX[s]]))];
    })));
  if (isHSC) {
    const e = hspc.engraftment(tr[tr.length - 1].y);
    console.log(`\nfold expansion     ${r2(e.fold_expansion, 1)}x`);
    console.log(`stem frequency     ${r2(e.stem_frequency, 3)}`);
    console.log(`homing capacity    ${r2(e.homing)}`);
    console.log(`ENGRAFTMENT        ${r2(e.engraftment, 2)}   (units x homing, per input cell)`);
  }
  const csv = arg('csv', null);
  if (csv && !isHSC) {
    console.log('\n->', writeTrajectory(csv, tr, (p) => M.observables(p.y)));
  }
} else {
  const { summary } = runPopulation(u, { n: cells, tEnd: T, h: 0.02, model: M, markers: MARKERS });
  console.log(table(['marker', '% positive', 'mean level', 'per input cell'],
    MARKERS.map((m) => [m, r2(summary[`pct_${m}`], 1), r2(summary[`mean_${m}`]),
      r2(summary[`yield_${m}`], 1) + ' %'])));
  console.log(`\nviability          ${r2(summary.viability)}`);
}
