// Canonical ordering of the state vector. Every module indexes through these
// names so the layout can change in one place.

export const SPECIES = [
  // --- Layer 0: drug ---
  'CHIR_in',   //  0  intracellular CHIR99021                 (uM)
  'V',         //  1  viability / fraction surviving          (0..1)

  // --- Layer 1: Wnt / beta-catenin ---
  'R',         //  2  activated Fzd/LRP6 fraction             (0..1)
  'AXIN2',     //  3  AXIN2 scaffold protein                  (a.u.)
  'Bcat_c',    //  4  free cytoplasmic beta-catenin           (nM)
  'Bcat_n',    //  5  free nuclear beta-catenin               (nM)
  'Bcat_p',    //  6  phospho-beta-catenin (degradation pool) (nM)
  'Bcat_m',    //  7  cadherin-bound beta-catenin             (nM)
  'TCF_B',     //  8  nuclear beta-catenin:TCF/LEF complex    (nM)  <- the Wnt OUTPUT
  'WNT_s',     //  9  autocrine WNT3/WNT5A in the medium       (nM)

  // --- Layer 2: fate GRN (all 0..1) ---
  'OCT4', 'SOX2', 'NANOG',
  'TBXT', 'EOMES', 'MIXL1', 'SNAI1',
  'MESP1', 'SOX17', 'GATA4', 'NKX25', 'TNNT2'
];

export const IDX = Object.fromEntries(SPECIES.map((s, i) => [s, i]));
export const N = SPECIES.length;

/** Undifferentiated hPSC (E8/mTeSR) at t = 0, no drug. */
export function initialState() {
  const y = new Float64Array(N);
  y[IDX.CHIR_in] = 0;
  y[IDX.V]       = 1;
  y[IDX.R]       = 0;
  y[IDX.AXIN2]   = 0.05;
  y[IDX.Bcat_c]  = 25;    // low free pool - destruction complex is active
  y[IDX.Bcat_n]  = 8;
  y[IDX.Bcat_p]  = 2;
  y[IDX.Bcat_m]  = 90;    // most beta-catenin is at the junctions in hPSC
  y[IDX.TCF_B]   = 1.0;
  y[IDX.WNT_s]   = 0;
  y[IDX.OCT4]    = 0.95;
  y[IDX.SOX2]    = 0.95;
  y[IDX.NANOG]   = 0.90;
  return y;
}

/** Convert a Float64Array state into a named object. */
export function toObj(y) {
  const o = {};
  for (let i = 0; i < N; i++) o[SPECIES[i]] = y[i];
  return o;
}
