import { P } from './params.js';
import { hill, rep, mm, pos } from './kinetics.js';
import { destructionComplex } from './interfaces.js';

// ---------------------------------------------------------------------------
// Layers 0-1, shared by every fate module.
//
// The first ten state variables are identical in every model built on this
// core, at identical indices, so this function can write into any of them.
// ---------------------------------------------------------------------------
export const SHARED_SPECIES = [
  'CHIR_in',   // 0  intracellular CHIR99021                 (uM)
  'V',         // 1  viability                                (0..1)
  'R',         // 2  activated Fzd/LRP6 fraction              (0..1)
  'AXIN2',     // 3  AXIN2 scaffold protein                   (a.u.)
  'Bcat_c',    // 4  free cytoplasmic beta-catenin            (nM)
  'Bcat_n',    // 5  free nuclear beta-catenin                (nM)
  'Bcat_p',    // 6  phospho-beta-catenin                     (nM)
  'Bcat_m',    // 7  cadherin-bound beta-catenin              (nM)
  'TCF_B',     // 8  nuclear beta-catenin:TCF/LEF             (nM)  <- the OUTPUT
  'WNT_s'      // 9  autocrine WNT ligand in the medium       (nM)
];
export const S = Object.fromEntries(SHARED_SPECIES.map((s, i) => [s, i]));

// ---------------------------------------------------------------------------
// CHIR99021 is ATP-competitive, so the dose a cell needs is set by
// Cheng-Prusoff against the cellular ATP pool, not by the enzymatic Ki:
//
//     Ki_app = Ki * (1 + [ATP]/Km_ATP) = 9.9 nM * (1 + 3000/15) ~ 2.0 uM
//
// a 201-fold shift, which is why a 6.7 nM inhibitor is pipetted at 3-12 uM.
// ---------------------------------------------------------------------------
export function gsk3Activity(chirIntracellular, par = P) {
  const c = par.chir;
  const KiApp = c.Ki_gsk3b * (1 + c.ATP_cell / c.Km_ATP);
  return { active: 1 / (1 + Math.pow(pos(chirIntracellular) / KiApp, c.hill)), KiApp };
}

export const cellularIC50 = (par = P) => gsk3Activity(0, par).KiApp / par.chir.Kp;

/**
 * Write d/dt for the shared layers into `d`.
 *
 * @param inp       medium composition from the protocol
 * @param fb        feedback from the fate module:
 *                    emt           0..1  junction dissolution (releases the cadherin pool)
 *                    wntSecretion  0..1  autocrine WNT ligand production
 * @returns derived quantities the fate module needs
 */
export function wntCore(y, inp, par, d, fb = { emt: 0, wntSecretion: 0 }) {
  const c = par.chir, w = par.wnt;
  const CHIR_in = y[S.CHIR_in], V = y[S.V], R = y[S.R], AXIN2 = y[S.AXIN2];
  const Bc = y[S.Bcat_c], Bn = y[S.Bcat_n], Bp = y[S.Bcat_p], Bm = y[S.Bcat_m];
  const W = y[S.TCF_B], WNT_s = y[S.WNT_s];

  // --- drug uptake and viability -------------------------------------------
  d[S.CHIR_in] = c.k_uptake * (c.Kp * inp.chir - CHIR_in);
  d[S.V] = -c.k_tox * hill(CHIR_in, c.IC50_tox, c.hill_tox) * V
           + c.k_recover * V * (1 - V);

  // --- receptor and autocrine ligand ---------------------------------------
  // IWP2 blocks palmitoylation so no NEW ligand is made; ligand already in the
  // medium still has to clear, which is why the inhibitor window is days long.
  d[S.WNT_s] = w.k_sec_wnt * fb.wntSecretion * (1 - inp.iwp) - w.k_clr_wnt * WNT_s;
  const occW = mm(inp.wnt3a + WNT_s, w.Kd_wnt);
  d[S.R] = w.k_on_R * occW * (1 - R) - w.k_off_R * R;

  // --- destruction complex -------------------------------------------------
  const aGSK3 = gsk3Activity(CHIR_in, par).active;
  const dcActivity = aGSK3 * (1 - w.eta_R * R);   // used for nuclear retention
  // AXIN2 is itself a Wnt target gene: the pathway's dominant negative feedback
  d[S.AXIN2] = w.k_syn_ax2 * hill(W, w.K_ax2, w.n_ax2) - w.k_deg_ax2 * AXIN2;
  // Explicit scaffold assembly, so a tankyrase inhibitor (IWR-1, which spares
  // Axin from degradation) is a different perturbation from a porcupine
  // inhibitor (IWP2, which blocks ligand). IWP2 cannot touch a CHIR response;
  // IWR-1 can, by raising total scaffold.
  const asm = destructionComplex(w.Axin0 + AXIN2, aGSK3, inp.iwr, par.dc);
  const DC = asm.active * (1 - w.eta_R * R);

  // --- beta-catenin pools ---------------------------------------------------
  const vPhos = w.k_phos * DC * mm(Bc, w.Km_b);
  // APC and Axin also shuttle beta-catenin OUT of the nucleus, so knocking the
  // complex down both spares it and retains it: two gains in series.
  const kExpEff = w.k_exp * (1 - w.exp_dc_gain + w.exp_dc_gain * dcActivity);
  const cdhTot = w.CDH1_0 * (1 - w.emt_loss * fb.emt);
  const vMem = w.k_mem_on * Bc * pos(cdhTot - Bm) - w.k_mem_off * Bm;
  const vTcf = w.k_tcf_on * Bn * pos(w.TCF_tot - W) - w.k_tcf_off * W;

  d[S.Bcat_c] = w.k_syn_b - vPhos - w.k_deg_b * Bc - w.k_imp * Bc + kExpEff * Bn - vMem;
  d[S.Bcat_n] = w.k_imp * Bc - kExpEff * Bn - vTcf;
  d[S.Bcat_p] = vPhos - w.k_deg_bp * Bp;
  d[S.Bcat_m] = vMem;
  d[S.TCF_B] = vTcf;

  return { aGSK3, dcActivity, DC, W, V };
}
