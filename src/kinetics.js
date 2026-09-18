// Shared rate-law helpers. x is clamped at >= 0 so a transient negative
// excursion from the integrator can never produce NaN from a fractional power.

export const pos = (x) => (x > 0 ? x : 0);

/** Activating Hill function: 0 -> 1 as x rises past K. */
export function hill(x, K, n) {
  const v = Math.pow(pos(x), n);
  return v / (Math.pow(K, n) + v);
}

/** Repressing Hill function: 1 -> 0 as x rises past K. */
export function rep(x, K, n) {
  const kn = Math.pow(K, n);
  return kn / (kn + Math.pow(pos(x), n));
}

/** Noisy-OR: either input alone can drive the output. */
export function orv(...xs) {
  let q = 1;
  for (const x of xs) q *= 1 - Math.min(1, pos(x));
  return 1 - q;
}

/** Michaelis-Menten saturation. */
export const mm = (x, Km) => pos(x) / (Km + pos(x));
