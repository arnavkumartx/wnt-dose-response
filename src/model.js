import { P } from './params.js';
import { IDX, N } from './state.js';
import { S, wntCore, gsk3Activity, cellularIC50 } from './wnt-core.js';
import { hill, rep, orv } from './kinetics.js';

export { gsk3Activity, cellularIC50 };

// ---------------------------------------------------------------------------
// Cardiac mesoderm module (hPSC -> cardiomyocyte).
//
// Layers 0-1 live in wnt-core.js and are shared with the cord blood module.
// Everything below is the fate layer.
// ---------------------------------------------------------------------------

const SCRATCH = new Float64Array(N);

export function rhs(t, y, u, par = P, out = SCRATCH) {
  const g = par.grn, w = par.wnt;
  const inp = u(t);
  const d = out; d.fill(0);

  const OCT4 = y[IDX.OCT4], SOX2 = y[IDX.SOX2], NANOG = y[IDX.NANOG];
  const TBXT = y[IDX.TBXT], EOMES = y[IDX.EOMES], MIXL1 = y[IDX.MIXL1], SNAI1 = y[IDX.SNAI1];
  const MESP1 = y[IDX.MESP1], SOX17 = y[IDX.SOX17], GATA4 = y[IDX.GATA4];
  const NKX25 = y[IDX.NKX25], TNNT2 = y[IDX.TNNT2];

  // Feedback from the fate layer back onto signalling:
  //   EMT dissolves junctions and releases the cadherin-bound beta-catenin pool;
  //   streak cells (TBXT) and the nascent mesoderm they become (GATA4) secrete
  //   the autocrine WNT that the inhibitor step exists to remove.
  const emt = hill(SNAI1, w.K_snai_cdh, w.n_snai_cdh);
  const wntSecretion = orv(hill(TBXT, w.K_endo_wnt, w.n_endo_wnt),
                           w.gata4_wnt * hill(GATA4, w.K_gata4_wnt, w.n_gata4_wnt));

  const { W } = wntCore(y, inp, par, d, { emt, wntSecretion });

  // --- fate GRN -------------------------------------------------------------
  const S_ = hill(inp.activin, g.K_act, g.n_act);   // SMAD2/3   (Activin A / Nodal)
  const B_ = hill(inp.bmp4,    g.K_bmp, g.n_bmp);   // SMAD1/5/8 (BMP4)

  const OS = hill(Math.sqrt(Math.max(OCT4, 0) * Math.max(SOX2, 0)), g.K_os, g.n_os);
  // Cells leave pluripotency through EITHER side of the streak. Keying this on
  // TBXT alone means SOX2 never falls in a high-Nodal endoderm run, where TBXT
  // stays low and EOMES/SOX17 carry the transition instead.
  const leaving = orv(TBXT, EOMES, SOX17);
  const exitPluri = rep(leaving, g.K_tbxt_rep, g.n_tbxt_rep);

  d[IDX.OCT4]  = g.k_oct4  * (g.leak + (1 - g.leak) * OS) - g.d_oct4 * OCT4;
  d[IDX.SOX2]  = g.k_sox2  * (g.leak + (1 - g.leak) * OS * exitPluri) - g.d_sox2 * SOX2;
  d[IDX.NANOG] = g.k_nanog * (g.leak + (1 - g.leak) * OS * exitPluri * orv(S_, g.nanog_basal))
                 - g.d_nanog * NANOG;

  // Cells that have adopted a definitive identity leave the streak state:
  // NKX2-5+ cardiac mesoderm and SOX17+ endoderm are both T-negative.
  const committed = orv(hill(NKX25, g.K_commit, g.n_commit), hill(SOX17, g.K_commit, g.n_commit));

  // TBXT: direct TCF/LEF target, also Nodal-inducible. SOX2 RAISES its Wnt
  // threshold rather than vetoing it - a hard veto makes the switch unflippable.
  d[IDX.TBXT] = g.k_tbxt * (g.leak + (1 - g.leak)
                * orv(hill(W, g.K_tcf_tbxt, g.n_tcf_tbxt), g.tbxt_act * S_)
                * (g.sox2_rep_min + (1 - g.sox2_rep_min) * rep(SOX2, g.K_sox2_rep, g.n_sox2_rep))
                * (1 - g.commit_rep * committed)) - g.d_tbxt * TBXT;

  // EOMES needs BOTH Wnt and strong SMAD2/3 - this AND-gate is the endoderm arm.
  d[IDX.EOMES] = g.k_eomes * (g.leak + (1 - g.leak) * S_ * hill(W, g.K_tcf_eomes, g.n_tcf_eomes))
                 - g.d_eomes * EOMES;

  d[IDX.MIXL1] = g.k_mixl1 * (g.leak + (1 - g.leak)
                 * hill(W, g.K_tcf_mixl, g.n_tcf_mixl) * orv(S_, B_)) - g.d_mixl1 * MIXL1;

  d[IDX.SNAI1] = g.k_snai1 * (g.leak + (1 - g.leak)
                 * orv(hill(W, g.K_tcf_snai, g.n_tcf_snai), TBXT)) - g.d_snai1 * SNAI1;

  const streak = hill(TBXT, g.K_str, g.n_str) * orv(hill(MIXL1, g.K_str, 2), 0.30);
  d[IDX.MESP1] = g.k_mesp1 * (g.leak + (1 - g.leak)
                 * streak * rep(SOX17, g.K_sox17_rep, g.n_sox17_rep)) - g.d_mesp1 * MESP1;

  // SOX17 and FOXA2 cross-activate into a self-reinforcing definitive-endoderm
  // network. The latch is what lets endoderm survive the CHIR washout: without
  // it SOX17 peaks on day 1 and collapses with the Wnt signal, and the model
  // makes no endoderm on the real protocol.
  d[IDX.SOX17] = g.k_sox17 * (g.leak + (1 - g.leak) * orv(
                   hill(EOMES, g.K_str, g.n_str) * S_,
                   g.sox17_self * hill(SOX17, g.K_sox17_self, g.n_sox17_self)))
                 - g.d_sox17 * SOX17;

  // GATA4 is the cardiac-mesoderm memory. MESP1 is a pulse that is gone within a
  // day; GATA4 latches by self-activation and carries the state into the Wnt-off
  // window. Without this latch the model cannot make cardiomyocytes on ANY
  // schedule, which is itself a testable prediction.
  d[IDX.GATA4] = g.k_gata4 * (g.leak + (1 - g.leak)
                 * orv(hill(MESP1, g.K_mesp1, g.n_mesp1),
                       g.gata4_self * hill(GATA4, g.K_gata4, g.n_gata4))) - g.d_gata4 * GATA4;

  // NKX2-5 requires cardiac mesoderm AND canonical Wnt to be OFF. This single
  // repression term is what forces the protocol to be biphasic. Its own
  // autoregulation is the commitment step: once a cell is a cardiac progenitor,
  // returning Wnt no longer reverses it.
  d[IDX.NKX25] = g.k_nkx25 * (g.leak + (1 - g.leak) * orv(
                   hill(GATA4, g.K_gata4, g.n_gata4) * rep(W, g.K_tcf_nkx, g.n_tcf_nkx),
                   g.nkx_self * hill(NKX25, g.K_nkx_self, g.n_nkx_self))) - g.d_nkx25 * NKX25;

  d[IDX.TNNT2] = g.k_tnnt2 * (g.leak + (1 - g.leak)
                 * hill(NKX25, g.K_nkx25, g.n_nkx25) * hill(GATA4, g.K_gata4, 2))
                 - g.d_tnnt2 * TNNT2;

  return d;
}

/** Derived read-outs that are useful but are not state variables. */
export function observables(y, par = P) {
  const aGSK3 = gsk3Activity(y[S.CHIR_in], par).active;
  return {
    GSK3_active: aGSK3,
    GSK3_inhib_pct: 100 * (1 - aGSK3),
    Bcat_total: y[S.Bcat_c] + y[S.Bcat_n] + y[S.Bcat_p] + y[S.Bcat_m] + y[S.TCF_B]
  };
}
