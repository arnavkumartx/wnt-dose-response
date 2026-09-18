// A protocol is a list of time windows, each holding the medium composition.
// Anything not named in a window falls back to BASAL.
//
//   chir    CHIR99021 in the medium                (uM)
//   wnt3a   recombinant WNT3A                      (nM)
//   iwp     porcupine inhibitor, 0..1              (IWP2 6 uM   -> ~0.95)
//   iwr     tankyrase inhibitor, 0..1              (IWR-1 5 uM  -> ~0.90)
//   activin Activin A / Nodal, normalised 0..1     (100 ng/mL  -> 1.0)
//   bmp4    BMP4, normalised 0..1                  (10 ng/mL   -> 0.6)

export const BASAL = { chir: 0, wnt3a: 0, iwp: 0, iwr: 0, activin: 0.18, bmp4: 0.10, cytokines: 1.0 };
// cytokines: SCF + TPO + FLT3L for the HSPC module, normalised 0..1.
// Only the cord blood model reads it; the cardiac model ignores it.
// ^ the non-zero activin/bmp4 floor is the autocrine NODAL/BMP tone that hPSC
//   colonies generate for themselves in E8/mTeSR, not added growth factor.

export function protocol(windows) {
  // Merged medium compositions are built once. The RHS is called ~10^5 times
  // per cell, so object spreads in here dominate the runtime if left inline.
  const sorted = [...windows]
    .sort((a, b) => a.from - b.from)
    .map((w) => ({ from: w.from, to: w.to, medium: Object.freeze({ ...BASAL, ...w.set }) }));
  const basal = Object.freeze({ ...BASAL });
  return (t) => {
    let cur = basal;
    for (let i = 0; i < sorted.length; i++) {
      const w = sorted[i];
      if (t >= w.from && t < w.to) cur = w.medium;
    }
    return cur;
  };
}

export const constant = (set) => protocol([{ from: -1e9, to: 1e9, set }]);

// --- Library of real protocols (t in hours from induction) -----------------

/** Maintenance: E8/mTeSR, no small molecules. */
export const P_MAINTENANCE = protocol([]);

/** CHIR pulse only: the standard mesoderm-induction arm. */
export const P_CHIR_PULSE = (doseUM, durationH) =>
  protocol([{ from: 0, to: durationH, set: { chir: doseUM } }]);

/** CHIR held for the whole run - the "why doesn't this work?" control. */
export const P_CHIR_CONTINUOUS = (doseUM) =>
  protocol([{ from: 0, to: 1e9, set: { chir: doseUM } }]);

/**
 * GiWi - Lian et al. 2012/2013: GSK3 inhibition then Wnt inhibition.
 * d0-d1  CHIR 6-12 uM  ->  primitive streak / cardiac mesoderm
 * d1-d3  washout       ->  MESP1 pulse resolves
 * d3-d5  IWP2/IWR-1    ->  canonical Wnt off, NKX2-5 / TNNT2 rise
 */
export const P_GIWI = ({ dose = 6, chirEnd = 24, iwpStart = 72, iwpEnd = 120 } = {}) =>
  protocol([
    { from: 0,        to: chirEnd, set: { chir: dose } },
    { from: iwpStart, to: iwpEnd,  set: { iwp: 0.95 } }
  ]);

/** Definitive endoderm (D'Amour / STEMdiff): CHIR spike + sustained Activin A. */
export const P_ENDODERM = ({ dose = 3, chirEnd = 24 } = {}) =>
  protocol([
    { from: 0, to: chirEnd, set: { chir: dose, activin: 1.0 } },
    { from: chirEnd, to: 72, set: { activin: 1.0 } }
  ]);
