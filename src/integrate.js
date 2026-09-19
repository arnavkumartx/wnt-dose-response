

/**
 * Classical RK4 with a fixed step. The fastest process in the model is
 * betaTrCP-mediated clearance of phospho-beta-catenin (6 /h, tau ~ 10 min),
 * so h = 0.005 h (18 s) sits far inside the stability region. 00_calibrate.js
 * verifies this by step halving.
 */
export function rk4(rhs, y0, { t0 = 0, tEnd = 168, h = 0.005, sample = 0.25, u, par }) {
  // Size is taken from the initial condition, not from any one model's
  // species list, so the same integrator drives the cardiac and HSPC models.
  const N = y0.length;
  let y = Float64Array.from(y0);
  let t = t0;
  const out = [{ t, y: Float64Array.from(y) }];
  let nextSample = t0 + sample;

  const k1 = new Float64Array(N), k2 = new Float64Array(N);
  const k3 = new Float64Array(N), k4 = new Float64Array(N);
  const tmp = new Float64Array(N);

  const steps = Math.round((tEnd - t0) / h);
  for (let s = 0; s < steps; s++) {
    let a = rhs(t, y, u, par);                       k1.set(a);
    for (let i = 0; i < N; i++) tmp[i] = y[i] + 0.5 * h * k1[i];
    a = rhs(t + 0.5 * h, tmp, u, par);               k2.set(a);
    for (let i = 0; i < N; i++) tmp[i] = y[i] + 0.5 * h * k2[i];
    a = rhs(t + 0.5 * h, tmp, u, par);               k3.set(a);
    for (let i = 0; i < N; i++) tmp[i] = y[i] + h * k3[i];
    a = rhs(t + h, tmp, u, par);                     k4.set(a);

    for (let i = 0; i < N; i++) {
      y[i] += (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
      if (!Number.isFinite(y[i])) {
        // A missing parameter reads as undefined and poisons one state silently:
        // pos() maps NaN to 0, so a Hill term quietly returns a valid-looking
        // number and the run completes with that regulator simply switched off.
        // Fail loudly instead.
        throw new Error('non-finite state at t=' + t.toFixed(3) + ' h, index ' + i
          + ' (check for a misspelled parameter name)');
      }
      if (y[i] < 0) y[i] = 0;                   // concentrations stay physical
    }
    t += h;

    if (t >= nextSample - 1e-9) {
      out.push({ t, y: Float64Array.from(y) });
      nextSample += sample;
    }
  }
  return out;
}

/** Integrate to steady state and return the final state vector. */
export function steadyState(rhs, y0, u, { tEnd = 400, h = 0.005, par } = {}) {
  const tr = rk4(rhs, y0, { tEnd, h, sample: tEnd, u, par });
  return tr[tr.length - 1].y;
}
