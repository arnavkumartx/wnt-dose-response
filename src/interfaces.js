import { pos } from './kinetics.js';

// ---------------------------------------------------------------------------
// Protein-protein interface layer.
//
// Binding equilibrates in seconds to minutes; transcription and fate take
// hours. So the interfaces are solved as a rapid equilibrium at each time
// point rather than carried as more ODEs - a standard quasi-steady-state
// reduction, not an approximation of convenience.
//
// Two interfaces matter here.
//
// 1. THE ARM GROOVE.  TCF/LEF, ICAT (CTNNBIP1), APC and E-cadherin all bind
//    OVERLAPPING surfaces on beta-catenin armadillo repeats ~3-9. They are
//    mutually exclusive, so they compete for one site, and the competition runs
//    through the shared pool of free beta-catenin rather than ligand-on-ligand.
//
// 2. THE BCL9/PYGO SITE.  BCL9 binds ARM1, a DIFFERENT surface, so it is not
//    competitive - it is required in addition. beta-catenin:TCF is therefore
//    not the same thing as transcriptional output: a complex without BCL9/Pygo
//    is assembled but quiet. This is why the model separates the two, and it is
//    what makes BCL9-site inhibitors mechanistically distinct from anything
//    that lowers beta-catenin levels.
// ---------------------------------------------------------------------------

/**
 * Competitive occupancy of one site on beta-catenin.
 *
 * Solves  Btot = Bfree + SUM_i Bfree*L_i/(Kd_i + Bfree)  for Bfree by bisection,
 * which is monotone in Bfree so bisection is unconditionally safe.
 *
 * @param Btot     total beta-catenin in the compartment (nM)
 * @param ligands  [{ name, total (nM), Kd (nM) }]
 */
export function armGrooveOccupancy(Btot, ligands) {
  const bound = (Bf) => ligands.reduce((a, L) => a + Bf * L.total / (L.Kd + Bf), 0);
  let lo = 0, hi = Math.max(Btot, 1e-9);
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if (mid + bound(mid) < Btot) lo = mid; else hi = mid;
  }
  const Bfree = 0.5 * (lo + hi);
  const out = { Bfree, complexes: {}, occupancy: {} };
  for (const L of ligands) {
    const c = Bfree * L.total / (L.Kd + Bfree);
    out.complexes[L.name] = c;
    out.occupancy[L.name] = L.total > 0 ? c / L.total : 0;   // fraction of ligand engaged
  }
  return out;
}

/** Defaults: no ICAT, BCL9 present. These reproduce the base model exactly. */
export const DEFAULT_INTERFACES = {
  TCF_tot: 15.0, Kd_TCF: 40.0,     // [M/A] TCF/LEF in the nucleus
  ICAT_tot: 0.0, Kd_ICAT: 12.0,    // [M]   ICAT binds ARM 10-12 + the TCF surface
  APCn_tot: 8.0, Kd_APCn: 120.0,   // [A]   nuclear APC, weak competitor + exporter
  BCL9_tot: 30.0, Kd_BCL9: 15.0,   // [A]   ARM1, NON-competitive co-activator
  bcl9_floor: 0.15                 // [F]   residual output with no BCL9 engaged
};

/**
 * Transcriptional output of the nuclear pool.
 *
 * Returns `W_eff`, the beta-catenin:TCF complex that is actually competent to
 * transcribe. With the defaults this is ~the raw complex, so the calibrated
 * cardiac and HSPC results are unchanged; it only departs when ICAT or a
 * BCL9-site inhibitor is introduced.
 */
export function nuclearOutput(Btot_nuclear, cfg = DEFAULT_INTERFACES) {
  const c = { ...DEFAULT_INTERFACES, ...cfg };
  const eq = armGrooveOccupancy(pos(Btot_nuclear), [
    { name: 'TCF', total: c.TCF_tot, Kd: c.Kd_TCF },
    { name: 'ICAT', total: c.ICAT_tot, Kd: c.Kd_ICAT },
    { name: 'APC', total: c.APCn_tot, Kd: c.Kd_APCn }
  ]);
  // BCL9 binds ARM1 on beta-catenin that is ALREADY in a TCF complex.
  const bcl9Frac = c.BCL9_tot / (c.Kd_BCL9 + c.BCL9_tot);
  const active = c.bcl9_floor + (1 - c.bcl9_floor) * bcl9Frac;
  return {
    ...eq,
    bcatTCF: eq.complexes.TCF,
    bcl9Fraction: bcl9Frac,
    W_eff: eq.complexes.TCF * active
  };
}

// ---------------------------------------------------------------------------
// Destruction complex assembly.
//
// Modelling the scaffold explicitly separates two inhibitors that the lumped
// version cannot tell apart:
//
//   IWP2  (porcupine)   blocks WNT palmitoylation -> no new LIGAND.
//                       Useless against CHIR, which acts below the receptor.
//   IWR-1 (tankyrase)   blocks AXIN1 PARsylation -> Axin is DEGRADED LESS, so
//                       total scaffold rises. That partially opposes CHIR,
//                       because more scaffold means more residual complex even
//                       when each GSK3 is mostly inhibited.
//
// That is a falsifiable difference: IWR-1 should blunt a CHIR response and
// IWP2 should not.
// ---------------------------------------------------------------------------
export const DEFAULT_DC = {
  Axin_tot: 1.0,      // a.u. normalised scaffold
  APC_tot: 100.0,     // nM  [M] APC is in large excess over Axin
  GSK3_tot: 50.0,     // nM  [M]
  CK1a_tot: 40.0,     // nM  [A]
  Kd_APC: 30.0,       // nM  [A] Axin:APC
  Kd_GSK3: 10.0,      // nM  [M] Axin:GSK3 is tight
  Kd_CK1: 25.0,       // nM  [A] Axin:CK1alpha
  tankyrase_gain: 2.5 // -   [F] fold Axin stabilisation at saturating IWR-1
};

/**
 * Fraction of Axin scaffold carrying a complete, catalytically competent
 * complex, and the total scaffold available.
 *
 * @param axinScaffold  Axin1 + AXIN2 in the same a.u. as the core model
 * @param gsk3Active    0..1 from the CHIR layer
 * @param iwr           0..1 tankyrase inhibitor occupancy
 */
export function destructionComplex(axinScaffold, gsk3Active, iwr = 0, cfg = DEFAULT_DC) {
  const c = { ...DEFAULT_DC, ...cfg };
  // Tankyrase inhibition spares Axin from degradation: more scaffold.
  const axin = axinScaffold * (1 + (c.tankyrase_gain - 1) * iwr);
  // APC, GSK3 and CK1alpha are all in excess over Axin, so each site fills
  // independently and the assembled fraction is the product of the three.
  const fAPC = c.APC_tot / (c.Kd_APC + c.APC_tot);
  // An ATP-competitive inhibitor does NOT stop GSK3 binding Axin - the enzyme
  // is still in the complex, just not turning over. So assembly uses TOTAL
  // GSK3, and the activity term multiplies the catalytic output instead.
  const fGSK = c.GSK3_tot / (c.Kd_GSK3 + c.GSK3_tot);
  const fCK1 = c.CK1a_tot / (c.Kd_CK1 + c.CK1a_tot);
  const assembled = fAPC * fGSK * fCK1;
  return { axin, fAPC, fGSK, fCK1, assembled, active: axin * assembled * gsk3Active };
}
