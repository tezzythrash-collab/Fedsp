#!/usr/bin/env python3
"""FeDSP fluid-demo oracle (Rank 4, Clay mapping).

Two solvers for the same 2D Taylor-Green / perturbed Navier-Stokes
problem (periodic [0, 2*pi)^2, vorticity-streamfunction, pseudo-spectral,
2/3 de-aliasing, RK4), run side by side:

  1. FIXED-POINT solver — mirrors the Fe code in ingots/fluid/src/lib.fe
     operation for operation, using the exact integer primitives from
     tools/gen_vectors.py (sat_mul, sat_add, sat_sub, sat_neg,
     q_reciprocal, from_ratio, fft_mirror, ifft_mirror). Its outputs
     are the exact expected values for the Fe tests: if Fe disagrees,
     the Fe code (or this mirror) is wrong.
  2. FLOAT solver — the reference implementation from
     ~/workspace/millennium/navier-stokes-2d/ (taylor_green_2d.run and
     perturbed_convergence.run, imported, not copied), in float64.
     The fixed-point run is measured against it: the round-off report.

Fixed-point conventions (identical to the Fe code; see lib.fe comments):
  - State is the spectral vorticity w_hat as a flat row-major [Cplx; N*N]
    array of Q15 raws, FFT ordering (k = 0..N/2-1, -N/2..-1).
  - Streamfunction: psi_hat = w_hat / K2 via Q15 reciprocal of K2
    (K2 = 0 mode forced to 0). Velocity: u_hat = i*ky*psi_hat,
    v_hat = -i*kx*psi_hat. Gradients: wx_hat = i*kx*w_hat,
    wy_hat = i*ky*w_hat. All complex-by-real products are sat_mul.
  - Nonlinear term: (u.wx + v.wy) in physical space, fft2'd, then the
    2/3 de-aliasing mask (zero modes with |kx| > N//3 or |ky| > N//3).
  - Diffusion: -nu*K2*w_hat with nu = 1/100, coefficient per mode
    from_ratio(K2, 100) in Q15.
  - Time step dt = 4/1000 (Q15 raw 131 = from_ratio(4, 1000)). RK4
    stage/full-step multiplies are done on raws with exact integer
    helpers (mul_dt, mul_dt_half, mul_dt_sixth) so the small constants
    dt/2 and dt/6 are represented exactly, not as coarse Q15 values.
  - Initial spectrum: analytic modes of the initial field, rounded to
    Q15 raws (float DFT of the analytic field on the grid, rounded to
    the nearest integer raw). TG: 4 modes; perturbed: 16 modes. The
    float solver starts from the analytic continuous field; the
    amplitude quantization (<= 2e-5 relative) is part of the measured
    difference.
  - Diagnostics per step: energy = 0.5*mean(u^2+v^2) from physical
    u, v (i64 saturating accumulation of sat_mul squares,
    round-half-up mean over the 2*nn squared terms — that mean IS
    the energy, the 0.5 is the division by 2*nn, no extra half);
    max|w| from physical w; BKM trapezoidal with dt/2 on raws.

Usage:
  python3 tools/fluid_oracle.py short    # only the 10-step N=32 run (fast)
  python3 tools/fluid_oracle.py full32   # + full 250-step N=32 runs
  python3 tools/fluid_oracle.py full64   # + full 250-step N=64 runs
  python3 tools/fluid_oracle.py all      # everything (slow in pure Python)
Writes tools/fluid_oracle.json (fixed-point trajectories + float
reference + round-off report).
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_vectors as g

import numpy as np

sys.path.insert(0, "/home/hatch/workspace/millennium/navier-stokes-2d")
import taylor_green_2d as tg
import perturbed_convergence as pc

I32_MAX = 2**31 - 1
I32_MIN = -(2**31)

DT_NUM = 4          # dt = 4/1000 = 0.004
DT_DEN = 1000
DT_RAW = g.from_ratio(DT_NUM, DT_DEN, 15)   # 131
NU_DEN = 100          # nu = 1/100 = 0.01


def clamp_i32(x):
    return max(I32_MIN, min(I32_MAX, x))


def trunc_div(a, b):
    q = abs(a) // abs(b)
    return q if (a < 0) == (b < 0) else -q


def mul_dt(x):
    """x * dt, dt = Q15 raw 131: (x*131 + 2**14) >> 15, saturating.
    Identical to sat_mul(x, 131, 15)."""
    return clamp_i32(((x * 131) + (1 << 14)) >> 15)


def mul_dt_half(x):
    """x * dt/2 = x*131/65536, round-half-up (ties toward +inf),
    saturating. Mirrors the Fe mul_dt_half helper exactly."""
    return clamp_i32(((x * 131) + (1 << 15)) >> 16)


def mul_dt_sixth(x):
    """x * dt/6 = x*131/196608, round-half-up (ties toward +inf),
    saturating. Mirrors the Fe mul_dt_sixth helper exactly."""
    num = x * 131
    half = 98304
    den = 196608
    if num >= 0:
        return clamp_i32(trunc_div(num + half, den))
    return clamp_i32(trunc_div(num - half, den))


def wavenumbers(n):
    """FFT-ordered wavenumbers k = 0..n/2-1, -n/2..-1."""
    return [i if i < n // 2 else i - n for i in range(n)]


def mode_tables(n):
    """Per-mode (flat row-major, idx = r*n + c, r -> ky, c -> kx):
    kx, ky, K2, inv_k2 (Q15 raw, 0 for the zero mode),
    diff coeff from_ratio(K2, 100) (Q15 raw), dealias mask."""
    k = wavenumbers(n)
    kx, ky, k2, inv, diff, deal = [], [], [], [], [], []
    kmax = n // 3
    for r in range(n):
        for c in range(n):
            x, y = k[c], k[r]
            k2v = x * x + y * y
            kx.append(x)
            ky.append(y)
            k2.append(k2v)
            inv.append(0 if k2v == 0 else g.q_reciprocal(k2v << 15, 15))
            diff.append(g.from_ratio(k2v, NU_DEN, 15))
            deal.append(abs(x) <= kmax and abs(y) <= kmax)
    return kx, ky, k2, inv, diff, deal


def fft2_mirror(a, n, tw_re, tw_im, rev):
    """2D FFT via row FFTs then column FFTs, mirroring the Fe fft2."""
    b = list(a)
    for r in range(n):
        row = b[r * n:(r + 1) * n]
        b[r * n:(r + 1) * n] = g.fft_mirror(row, tw_re, tw_im, rev)
    for c in range(n):
        col = [b[r * n + c] for r in range(n)]
        col = g.fft_mirror(col, tw_re, tw_im, rev)
        for r in range(n):
            b[r * n + c] = col[r]
    return b


def ifft2_mirror(a, n, tw_re, tw_im, rev):
    b = list(a)
    for r in range(n):
        row = b[r * n:(r + 1) * n]
        b[r * n:(r + 1) * n] = g.ifft_mirror(row, tw_re, tw_im, rev)
    for c in range(n):
        col = [b[r * n + c] for r in range(n)]
        col = g.ifft_mirror(col, tw_re, tw_im, rev)
        for r in range(n):
            b[r * n + c] = col[r]
    return b


def cplx_mul_real(a, s):
    """Complex Q15 (re, im) times a real Q15 raw s: both sat_mul."""
    return (g.sat_mul(a[0], s, 15), g.sat_mul(a[1], s, 15))


def cplx_mul_i(a, k_raw):
    """Multiply complex a by i*k where k is a Q15 raw holding an
    integer k: (re, im) -> (-im*k, re*k). Negation is saturating
    (sat_neg), products are sat_mul — mirrors the Fe code."""
    return (g.sat_mul(g.sat_neg(a[1]), k_raw, 15),
            g.sat_mul(a[0], k_raw, 15))


def cplx_mul_neg_i(a, k_raw):
    """Multiply complex a by -i*k: (re, im) -> (im*k, -re*k)."""
    return (g.sat_mul(a[1], k_raw, 15),
            g.sat_mul(g.sat_neg(a[0]), k_raw, 15))


def cplx_neg(a):
    """Componentwise saturating negation."""
    return (g.sat_neg(a[0]), g.sat_neg(a[1]))


def init_spectrum(field, n):
    """Analytic initial spectrum as Q15 raws: float DFT (numpy) of the
    analytic field on the grid, rounded to the nearest integer raw.
    Sparse: only the analytic modes survive rounding (TG: 4 modes,
    perturbed: 16 modes); everything else rounds to exactly 0."""
    xs = []
    for r in range(n):
        for c in range(n):
            x = 2.0 * math.pi * c / n
            y = 2.0 * math.pi * r / n
            if field == "tg":
                w = 2.0 * math.sin(x) * math.sin(y)
            else:
                w = (2.0 * math.sin(x) * math.sin(y)
                     + 0.5 * math.sin(2 * x) * math.cos(3 * y)
                     + 0.3 * math.cos(3 * x) * math.sin(y)
                     - 0.2 * math.sin(x) * math.cos(2 * y))
            xs.append(w)
    arr = np.array(xs, dtype=np.float64).reshape(n, n)
    spec = np.fft.fft2(arr)
    out = []
    for v in spec.flat:
        # DFT coefficients are in field units; Q15 raws scale by 32768.
        re, im = int(round(v.real * 32768)), int(round(v.imag * 32768))
        out.append((re, im))
    return out


class FixedSolver:
    """Fixed-point 2D spectral solver, mirroring the Fe code exactly."""

    def __init__(self, n, field):
        self.n = n
        self.nn = n * n
        self.tw_re, self.tw_im = g.twiddles(n)
        self.rev = g.bitrev(n)
        self.kx, self.ky, self.k2, self.inv_k2, self.diff, self.deal = mode_tables(n)
        self.kx_raw = [x << 15 for x in self.kx]
        self.ky_raw = [y << 15 for y in self.ky]
        self.w_hat = init_spectrum(field, n)
        self.bkm = 0  # Q15 raw
        self.diag = []  # (energy_raw, wmax_raw, bkm_raw) per step incl. t=0

    def rhs(self, wh):
        n, nn = self.n, self.nn
        psi = [cplx_mul_real(wh[i], self.inv_k2[i]) for i in range(nn)]
        u_hat = [cplx_mul_i(psi[i], self.ky_raw[i]) for i in range(nn)]
        v_hat = [cplx_mul_neg_i(psi[i], self.kx_raw[i]) for i in range(nn)]
        wx_hat = [cplx_mul_i(wh[i], self.kx_raw[i]) for i in range(nn)]
        wy_hat = [cplx_mul_i(wh[i], self.ky_raw[i]) for i in range(nn)]
        u = [p[0] for p in ifft2_mirror(u_hat, n, self.tw_re, self.tw_im, self.rev)]
        v = [p[0] for p in ifft2_mirror(v_hat, n, self.tw_re, self.tw_im, self.rev)]
        wx = [p[0] for p in ifft2_mirror(wx_hat, n, self.tw_re, self.tw_im, self.rev)]
        wy = [p[0] for p in ifft2_mirror(wy_hat, n, self.tw_re, self.tw_im, self.rev)]
        nl = [(g.sat_add_i32(g.sat_mul(u[i], wx[i], 15),
                              g.sat_mul(v[i], wy[i], 15)), 0) for i in range(nn)]
        nl_hat = fft2_mirror(nl, n, self.tw_re, self.tw_im, self.rev)
        out = []
        for i in range(nn):
            if not self.deal[i]:
                out.append((0, 0))
                continue
            d = cplx_mul_real(wh[i], self.diff[i])
            neg = cplx_neg(nl_hat[i])
            out.append((g.sat_sub_i32(neg[0], d[0]),
                        g.sat_sub_i32(neg[1], d[1])))
        return out

    def add_scaled(self, wh, k, f):
        """wh + f(k) componentwise, f in {mul_dt, mul_dt_half}."""
        return [(g.sat_add_i32(wh[i][0], f(k[i][0])),
                 g.sat_add_i32(wh[i][1], f(k[i][1]))) for i in range(self.nn)]

    def step(self):
        nn = self.nn
        w = self.w_hat
        k1 = self.rhs(w)
        k2 = self.rhs(self.add_scaled(w, k1, mul_dt_half))
        k3 = self.rhs(self.add_scaled(w, k2, mul_dt_half))
        k4 = self.rhs(self.add_scaled(w, k3, mul_dt))
        s = []
        for i in range(nn):
            # S = k1 + 2*k2 + 2*k3 + k4, doubling via saturating add
            # (mirrors the Fe code, which doubles with sat_add).
            s.append((g.sat_add_i32(g.sat_add_i32(k1[i][0], g.sat_add_i32(k2[i][0], k2[i][0])),
                                     g.sat_add_i32(g.sat_add_i32(k3[i][0], k3[i][0]), k4[i][0])),
                      g.sat_add_i32(g.sat_add_i32(k1[i][1], g.sat_add_i32(k2[i][1], k2[i][1])),
                                     g.sat_add_i32(g.sat_add_i32(k3[i][1], k3[i][1]), k4[i][1]))))
        self.w_hat = [(g.sat_add_i32(w[i][0], mul_dt_sixth(s[i][0])),
                       g.sat_add_i32(w[i][1], mul_dt_sixth(s[i][1])))
                      for i in range(nn)]

    def diagnostics(self):
        n, nn = self.n, self.nn
        psi = [cplx_mul_real(self.w_hat[i], self.inv_k2[i]) for i in range(nn)]
        u_hat = [cplx_mul_i(psi[i], self.ky_raw[i]) for i in range(nn)]
        v_hat = [cplx_mul_neg_i(psi[i], self.kx_raw[i]) for i in range(nn)]
        u = [p[0] for p in ifft2_mirror(u_hat, n, self.tw_re, self.tw_im, self.rev)]
        v = [p[0] for p in ifft2_mirror(v_hat, n, self.tw_re, self.tw_im, self.rev)]
        w = [p[0] for p in ifft2_mirror(self.w_hat, n, self.tw_re, self.tw_im, self.rev)]
        # energy = 0.5 * mean(u^2 + v^2): Q15 products (sat_mul, like
        # dot) accumulated in i64 with saturation, round-half-up mean
        # over the 2*nn squared terms. That mean IS the energy: the
        # 0.5 is the division by 2*nn rather than nn (no extra half).
        acc = 0
        for i in range(nn):
            acc = g.clamp_i64(acc + g.sat_mul(u[i], u[i], 15))
            acc = g.clamp_i64(acc + g.sat_mul(v[i], v[i], 15))
        m = 2 * nn
        q = acc // m
        rem = acc - q * m
        if rem >= (m + 1) // 2:
            q += 1
        energy_raw = g.clamp_i32(q)
        wmax_raw = max(abs(x) for x in w)
        return energy_raw, wmax_raw

    def run(self, steps):
        e0, w0 = self.diagnostics()
        self.diag.append((e0, w0, self.bkm))
        w_prev = w0
        for _ in range(steps):
            self.step()
            e, wmax = self.diagnostics()
            s = g.sat_add_i32(w_prev, wmax)
            self.bkm = g.sat_add_i32(self.bkm, mul_dt_half(s))
            self.diag.append((e, wmax, self.bkm))
            w_prev = wmax
        return self.diag


def q15_to_float(raw):
    return raw / 32768.0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "short"
    configs = [("tg", 32, 10), ("perturbed", 32, 10)]
    if mode in ("full32", "all"):
        configs += [("tg", 32, 250), ("perturbed", 32, 250)]
    if mode in ("full64", "all"):
        configs += [("tg", 64, 250), ("perturbed", 64, 250)]
    out = {"dt_raw": DT_RAW, "dt": DT_NUM / DT_DEN, "nu": 1.0 / NU_DEN}
    for field, n, steps in configs:
        key = f"{field}_n{n}_{steps}steps"
        print(f"--- {key} ---", flush=True)
        fs = FixedSolver(n, field)
        diag = fs.run(steps)
        e_f, w_f = diag[-1][0], diag[-1][1]
        bkm_f = diag[-1][2]
        mono = all(diag[i + 1][0] <= diag[i][0] for i in range(len(diag) - 1))
        t_end = steps * DT_NUM / DT_DEN
        if field == "tg":
            fr = tg.run(N=n, nu=0.01, t_end=t_end, dt=DT_NUM / DT_DEN)
            e_ref, w_ref, bkm_ref = fr["energy_final"], fr["wmax_final"], fr["bkm_integral"]
        else:
            fr = pc.run(N=n, nu=0.01, t_end=t_end, dt=DT_NUM / DT_DEN)
            e_ref, w_ref, bkm_ref = fr["energy_final"], fr["wmax_final"], fr["bkm_integral"]
        e_fix, w_fix = q15_to_float(e_f), q15_to_float(w_f)
        bkm_fix = q15_to_float(bkm_f)
        out[key] = {
            "field": field, "N": n, "steps": steps, "t_end": t_end,
            "fixed": {"energy_raw": e_f, "energy": e_fix,
                      "wmax_raw": w_f, "wmax": w_fix,
                      "bkm_raw": bkm_f, "bkm": bkm_fix,
                      "energy_monotone": mono,
                      "diag": diag},
            "float": {"energy": e_ref, "wmax": w_ref, "bkm": bkm_ref},
            "roundoff": {
                "energy_rel": abs(e_fix - e_ref) / e_ref,
                "wmax_rel": abs(w_fix - w_ref) / w_ref,
                "bkm_rel": abs(bkm_fix - bkm_ref) / bkm_ref,
            },
        }
        print(f"  fixed E={e_fix:.8f} wmax={w_fix:.6f} BKM={bkm_fix:.4f} mono={mono}")
        print(f"  float E={e_ref:.8f} wmax={w_ref:.6f} BKM={bkm_ref:.4f}")
        r = out[key]["roundoff"]
        print(f"  rel err: E={r['energy_rel']:.3e} wmax={r['wmax_rel']:.3e} BKM={r['bkm_rel']:.3e}")
    with open("tools/fluid_oracle.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote tools/fluid_oracle.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
