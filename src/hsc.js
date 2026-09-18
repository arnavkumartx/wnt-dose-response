import { P } from './params.js';
import { SHARED_SPECIES, S, wntCore, gsk3Activity } from './wnt-core.js';
import { hill, rep, orv } from './kinetics.js';

// ---------------------------------------------------------------------------
// Cord blood CD34+ HSPC module.
//
// The point of this module is that in haematopoiesis the Wnt dose-response is
// NOT monotonic, and the model is built so that falls out rather than being
// asserted. Two direct Wnt target genes have different thresholds:
//
//     HOXB4  switches on at LOW  beta-catenin  -> self-renewal
//     MYC    switches on at HIGH beta-catenin  -> niche exit, differentiation
//
// so mild pathway activation raises self-renewal while strong activation
// overruns it. That is the mechanism behind the "dosage-dependent" picture in
// Luis 2011, and it is why gain-of-function studies report HSC expansion or HSC
// depletion depending on how hard they pushed.
//
// Note what this implies for CHIR: in cord blood HSPC, GSK3 inhibition is a
// self-renewal / engraftment intervention, not a differentiation one. The
// published cord blood result is *delayed* expansion with *preserved* stem cell
// activity (Holmes 2008), and the same compound given in vivo after transplant
// makes regeneration worse (Shen 2014).
// ---------------------------------------------------------------------------

export const SPECIES = [
  ...SHARED_SPECIES,
  'HOXB4',   // 10  self-renewal arm, low-threshold Wnt target
  'MYC',     // 11  differentiation arm, high-threshold Wnt target
  'CDKN1C',  // 12  p57 - quiescence, up with GSK3 inhibition
  'CCND1',   // 13  cyclin D1 - cycling, down with GSK3 inhibition
  'CXCR4',   // 14  stromal adherence / marrow homing
  'GATA2',   // 15  HSC maintenance
  'GATA1',   // 16  erythro-megakaryocytic output
  'SPI1',    // 17  PU.1, myeloid output
  'LYMPH',   // 18  lymphoid competence (blocked by high Wnt)
  'CD34', 'CD90', 'CD45RA',  // 19-21  immunophenotype
  'STEM',    // 22  repopulating potential per cell (the SRC unit)
  'Ncell'    // 23  fold expansion relative to input
];
export const IDX = Object.fromEntries(SPECIES.map((s, i) => [s, i]));
export const N = SPECIES.length;

/** Freshly thawed, CD34-selected cord blood cells at t = 0. */
export function initialState() {
  const y = new Float64Array(N);
  y[S.V] = 1; y[S.AXIN2] = 0.05;
  y[S.Bcat_c] = 25; y[S.Bcat_n] = 4.3; y[S.Bcat_p] = 2;
  y[S.Bcat_m] = 95; y[S.TCF_B] = 1.2;
  y[IDX.HOXB4] = 0.20; y[IDX.MYC] = 0.10; y[IDX.GATA2] = 0.80;
  y[IDX.CDKN1C] = 0.60; y[IDX.CCND1] = 0.15;   // quiescent on arrival
  y[IDX.CXCR4] = 0.55;
  y[IDX.CD34] = 0.95; y[IDX.CD90] = 0.55; y[IDX.CD45RA] = 0.10;
  y[IDX.STEM] = 1.00;                           // full repopulating potential
  y[IDX.Ncell] = 1.00;                          // fold expansion
  return y;
}

const SCRATCH = new Float64Array(N);

export function rhs(t, y, u, par = P, out = SCRATCH) {
  const h = par.hsc;
  const inp = u(t);
  const d = out; d.fill(0);

  const HOXB4 = y[IDX.HOXB4], MYC = y[IDX.MYC], p57 = y[IDX.CDKN1C];
  const CCND1 = y[IDX.CCND1], CXCR4 = y[IDX.CXCR4], GATA2 = y[IDX.GATA2];
  const GATA1 = y[IDX.GATA1], SPI1 = y[IDX.SPI1], LYMPH = y[IDX.LYMPH];
  const CD34 = y[IDX.CD34], CD90 = y[IDX.CD90], CD45RA = y[IDX.CD45RA];
  const STEM = y[IDX.STEM], Ncell = y[IDX.Ncell];

  // HSPC in suspension keep their junction-free state and do not run the
  // autocrine WNT loop that primitive-streak cells do.
  const { W, V } = wntCore(y, inp, par, d, { emt: 0, wntSecretion: 0 });

  const cyt = hill(inp.cytokines, h.K_cyt, h.n_cyt);   // SCF + TPO + FLT3L

  // --- two direct Wnt targets, different thresholds -------------------------
  d[IDX.HOXB4] = h.k_hoxb4 * (h.leak + (1 - h.leak) * hill(W, h.K_w_hoxb4, h.n_w_hoxb4))
                 - h.d_hoxb4 * HOXB4;
  d[IDX.MYC]   = h.k_myc * (h.leak + (1 - h.leak)
                 * orv(hill(W, h.K_w_myc, h.n_w_myc), 0.45 * cyt)) - h.d_myc * MYC;

  // --- cell cycle: GSK3 inhibition lengthens cycle time (Ko 2011) -----------
  d[IDX.CDKN1C] = h.k_cdkn1c * (h.leak + (1 - h.leak) * hill(W, h.K_w_p57, h.n_w_p57))
                  - h.d_cdkn1c * p57;
  d[IDX.CCND1]  = h.k_ccnd1 * (h.leak + (1 - h.leak) * cyt * rep(W, h.K_w_ccnd1, h.n_w_ccnd1))
                  - h.d_ccnd1 * CCND1;
  d[IDX.CXCR4]  = h.k_cxcr4 * (h.leak + (1 - h.leak) * hill(W, h.K_w_cxcr4, h.n_w_cxcr4))
                  - h.d_cxcr4 * CXCR4;

  const prolif = h.k_prolif * CCND1 * rep(p57, h.K_p57, h.n_p57) * V;
  d[IDX.Ncell] = prolif * Ncell;

  // --- self-renewal vs differentiation --------------------------------------
  const selfRenew = hill(HOXB4, h.K_hoxb4, h.n_hoxb4) * rep(MYC, h.K_myc_rep, h.n_myc_rep);
  const diffDrive = hill(MYC, h.K_diff, h.n_diff) * cyt;

  d[IDX.GATA2] = h.k_gata2 * (h.leak + (1 - h.leak) * selfRenew) - h.d_gata2 * GATA2;

  // Repopulating potential relaxes toward the self-renewal set point, and it
  // relaxes FASTER the faster the cells cycle: every division is a chance to
  // differentiate. This is why unmanipulated ex vivo expansion trades cell
  // number for engraftment.
  const dStemEff = h.d_stem * (1 + h.k_div_loss * prolif / h.k_prolif);
  d[IDX.STEM] = dStemEff * (selfRenew - STEM);

  // --- lineage output --------------------------------------------------------
  d[IDX.GATA1] = h.k_gata1 * (h.leak + (1 - h.leak) * diffDrive * rep(SPI1, h.K_ant, h.n_ant))
                 - h.d_gata1 * GATA1;
  d[IDX.SPI1]  = h.k_spi1 * (h.leak + (1 - h.leak) * diffDrive * rep(GATA1, h.K_ant, h.n_ant))
                 - h.d_spi1 * SPI1;
  // High canonical Wnt suppresses lymphoid development (Shen 2014: T-cell
  // depletion persisted even after the inhibitor was withdrawn).
  d[IDX.LYMPH] = h.k_lymph * (h.leak + (1 - h.leak) * diffDrive * rep(W, h.K_w_lymph, h.n_w_lymph))
                 - h.d_lymph * LYMPH;

  // --- immunophenotype -------------------------------------------------------
  const committed = orv(hill(GATA1, h.K_diff, h.n_diff), hill(SPI1, h.K_diff, h.n_diff));
  d[IDX.CD34]   = h.k_cd34 * (h.leak + (1 - h.leak)
                  * orv(hill(STEM, h.K_stem, h.n_stem), hill(GATA2, h.K_diff, 2))
                  * (1 - 0.8 * committed)) - h.d_cd34 * CD34;
  d[IDX.CD90]   = h.k_cd90 * (h.leak + (1 - h.leak) * hill(STEM, h.K_stem, h.n_stem))
                  - h.d_cd90 * CD90;
  d[IDX.CD45RA] = h.k_cd45ra * (h.leak + (1 - h.leak) * committed) - h.d_cd45ra * CD45RA;

  return d;
}

/**
 * Engraftment read-outs. `units` is repopulating units harvested per input
 * cell; `engraftment` additionally weights each unit by its homing capacity,
 * which is the axis Ko 2011 found GSK3 inhibition acting on - more output per
 * unit rather than more units.
 */
export function engraftment(y, par = P) {
  const h = par.hsc;
  const homing = h.homing_floor + (1 - h.homing_floor) * hill(y[IDX.CXCR4], h.K_cxcr4, h.n_cxcr4);
  const units = y[IDX.Ncell] * y[IDX.STEM] * y[S.V];
  return {
    fold_expansion: y[IDX.Ncell],
    stem_frequency: y[IDX.STEM],
    homing,
    units,
    engraftment: units * homing,
    GSK3_active: gsk3Activity(y[S.CHIR_in], par).active
  };
}
