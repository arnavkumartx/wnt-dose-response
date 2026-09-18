import fs from 'node:fs';
import path from 'node:path';
import { SPECIES } from './state.js';

const OUT = path.resolve(process.cwd(), 'out');

export function writeCSV(name, rows, header) {
  fs.mkdirSync(OUT, { recursive: true });
  const f = path.join(OUT, name);
  const lines = [header.join(','), ...rows.map((r) => r.map(fmt).join(','))];
  fs.writeFileSync(f, lines.join('\n') + '\n');
  return f;
}

export function writeTrajectory(name, traj, extra = () => ({})) {
  const extraKeys = Object.keys(extra(traj[0]));
  const header = ['t_h', ...SPECIES, ...extraKeys];
  const rows = traj.map((p) => {
    const e = extra(p);
    return [p.t, ...Array.from(p.y), ...extraKeys.map((k) => e[k])];
  });
  return writeCSV(name, rows, header);
}

const fmt = (v) => (typeof v === 'number' ? (Number.isFinite(v) ? v.toPrecision(6) : '') : v);

/** Fixed-width table for terminal output. */
export function table(header, rows) {
  const all = [header, ...rows.map((r) => r.map(String))];
  const w = header.map((_, i) => Math.max(...all.map((r) => (r[i] ?? '').length)));
  const line = (r) => r.map((c, i) => String(c ?? '').padStart(w[i])).join('  ');
  return [line(header), w.map((n) => '-'.repeat(n)).join('  '), ...rows.map((r) => line(r.map(String)))].join('\n');
}

export const r2 = (x, n = 2) => Number(x).toFixed(n);
