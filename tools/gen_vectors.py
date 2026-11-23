#!/usr/bin/env python3
"""FeDSP test-vector oracle.

Computes every expected value embedded in the Fe tests using exact
arbitrary-precision integer arithmetic (Python ints never overflow), so
the Fe test vectors are independent of the Fe implementation itself.

Semantics mirrored exactly from src/lib.fe (documented there too):
  - Q<FRAC> value = raw / 2**FRAC, raw is a signed 32-bit integer.
  - sat_mul: (a*b + 2**(FRAC-1)) >> FRAC, then clamp to i32 range.
    (Round-half-up: ties round toward +inf; matches Fe arithmetic >>.)
  - sat_div: trunc((a << FRAC) / b) toward zero, clamp to i32;
    divisor 0 -> MAX (a>0), MIN (a<0), 0 (a==0).
  - from_int / from_i64: saturating, via i64 range check before shifting.
  - from_ratio(num, den): (num << FRAC) / den with round-to-nearest
    (ties away from zero), den == 0 -> MAX/MIN/0 by sign of num.
    The rounding adjustment follows the sign of the numerator, so
    negative denominators round correctly (v0.4.1 audit fix).
  - Hann window: round(0.5 * (1 - cos(2*pi*n/(N-1))) * 32767), Q15 raw.
  - dot (v0.2.0): i64 accumulation with saturation at every step
    (clamp_i64 after each add), then round-to-nearest and clamp to i32.
  - convert_q: Q<FROM> -> Q<TO> by shift, rounding to nearest
    (ties toward +inf) when shifting right, saturating.
  - biquad: transposed direct form II, one sat_mul/sat_add/sat_sub per
    term, in the same order as Biquad::step.
  - fir: delay-line shift, then the saturating dot.
  - mean: saturating i64 sum, round-to-nearest division by N
    (ties toward +inf). rms: saturating sum of squares, integer sqrt
    of the mean, rounded to nearest. energy: dot(x, x).
  - q_sqrt: isqrt(raw << FRAC), rounded to nearest; negative -> 0.
  - q_reciprocal: exactly rounded 1/x by one i64 division
    (2**(2*FRAC) / raw, ties toward +inf), saturating; 0 ulp on
    every input (v0.4.1: replaced the 5-step Newton iteration, which
    was ~1 ulp of the normalized reciprocal but up to 2**k ulp
    after rescaling).
  - differential (v0.4.0): seeded random inputs (seed 20261120,
    consumed in the fixed order dot -> fir -> biquad -> fft -> scalars
    -> demo-contract frame) through dot, FIR, biquad, FFT, sqrt,
    reciprocal, and convert_q. Differential coverage on top of the
    fixed vectors above, not a replacement: every expected value is
    computed here, never copied from Fe output.

Outputs:
  - tools/vectors.json : every vector used by the Fe tests.
  - tools/hann_tables.txt : Fe literals for HANN16/32/64.
  - tools/window_tables.txt : Fe literals for HAMMING16/32/64 and
    BLACKMAN16/32/64.
Run: python3 tools/gen_vectors.py
"""
import json
import math
import random

I32_MAX = 2**31 - 1
I32_MIN = -(2**31)
I64_MAX = 2**63 - 1
I64_MIN = -(2**63)


def clamp_i32(x):
    return max(I32_MIN, min(I32_MAX, x))


def clamp_i64(x):
    return max(I64_MIN, min(I64_MAX, x))


def sat_add_i32(a, b):
    return clamp_i32(a + b)


def sat_sub_i32(a, b):
    return clamp_i32(a - b)


def isqrt_round(x):
    """Integer square root, rounded to nearest (ties toward +inf).

    Mirrors the Fe Newton iteration in src/lib.fe: floor sqrt, then round
    up iff the remainder exceeds the root (r > s ⟺ M is closer to (s+1)^2).
    """
    if x <= 0:
        return 0
    s = math.isqrt(x)
    r = x - s * s
    return s + 1 if r > s else s


def trunc_div(a, b):
    """Integer division truncating toward zero (Fe i64 `/` semantics)."""
    q = abs(a) // abs(b)
    return q if (a < 0) == (b < 0) else -q


def sat_mul(a, b, frac):
    v = a * b
    if frac >= 1:
        v = (v + (1 << (frac - 1))) >> frac  # Python >> is arithmetic
    return clamp_i32(v)


def sat_div(a, b, frac):
    if b == 0:
        return I32_MAX if a > 0 else (I32_MIN if a < 0 else 0)
    return clamp_i32(trunc_div(a << frac, b))


def convert_q(raw, frm, to):
    """Q<FROM> raw -> Q<TO> raw, rounding to nearest (ties toward +inf),
    saturating. Mirrors convert_q in src/lib.fe exactly."""
    if frm == to:
        return clamp_i32(raw)
    if frm > to:
        d = frm - to
        return clamp_i32((raw + (1 << (d - 1))) >> d)
    return clamp_i32(raw << (to - frm))


def q_sqrt(raw, frac):
    """Q<FRAC> square root, raw in -> raw out, via integer sqrt of
    raw << FRAC, rounded to nearest. Negative input -> 0.
    Mirrors Q::sqrt in src/lib.fe."""
    if raw <= 0:
        return 0
    return clamp_i32(isqrt_round(raw << frac))


def q_reciprocal(raw, frac):
    """Q<FRAC> reciprocal, raw in -> raw out, exactly rounded to
    nearest (ties toward +inf), saturating. One i64 division:
    1/x in Q<FRAC> has raw 2^(2*FRAC) / raw. Mirrors Q::reciprocal
    in src/lib.fe exactly (v0.4.1: the old 5-step Newton mirror is
    gone — it was ~1 ulp of the normalized reciprocal but up to
    2^k ulp after rescaling)."""
    if raw == 0:
        return I32_MAX
    if frac == 0:
        return raw if raw in (1, -1) else 0
    if frac >= 31:
        return I32_MIN if raw < 0 else I32_MAX
    num = 1 << (2 * frac)
    den = raw
    if den > 0:
        return clamp_i32(trunc_div(num + den // 2, den))
    m = -den
    q, r = num // m, num % m
    return clamp_i32(-(q + 1) if 2 * r > m else -q)


def from_int(x, frac):
    max_int = I32_MAX >> frac
    min_int = I32_MIN >> frac  # arithmetic shift: floor, matches Fe
    if x > max_int:
        return I32_MAX
    if x < min_int:
        return I32_MIN
    return x << frac


def from_ratio(num, den, frac):
    if den == 0:
        return I32_MAX if num > 0 else (I32_MIN if num < 0 else 0)
    v = num << frac
    # round to nearest, ties away from zero. The adjustment follows
    # the sign of the *numerator* (v0.4.1 audit fix: the old branch
    # keyed on whether v and den agreed in sign, which rounded
    # toward zero for every negative denominator).
    half = abs(den) // 2
    if v >= 0:
        return clamp_i32(trunc_div(v + half, den))
    else:
        return clamp_i32(trunc_div(v - half, den))


def hann(n_points):
    return [
        int(round(0.5 * (1.0 - math.cos(2.0 * math.pi * n / (n_points - 1))) * 32767))
        for n in range(n_points)
    ]


def hamming(n_points):
    return [
        int(round((0.54 - 0.46 * math.cos(2.0 * math.pi * n / (n_points - 1))) * 32767))
        for n in range(n_points)
    ]


def blackman(n_points):
    return [
        int(round((0.42 - 0.5 * math.cos(2.0 * math.pi * n / (n_points - 1))
                   + 0.08 * math.cos(4.0 * math.pi * n / (n_points - 1))) * 32767))
        for n in range(n_points)
    ]


def dot(a, b, frac):
    """Saturating dot product: i64 accumulation with saturation at every
    step, then round-to-nearest (ties toward +inf) and clamp to i32.
    Mirrors dot in src/lib.fe exactly (v0.2.0 fix: the accumulator used
    to be a plain i64 add, which reverted on overflow)."""
    acc = 0
    for x, y in zip(a, b):
        acc = clamp_i64(acc + x * y)
    if frac >= 1:
        acc = clamp_i64(acc + (1 << (frac - 1))) >> frac
    return clamp_i32(acc)


def mean_py(vals):
    """Mean of Q15 raws, rounding to nearest (ties toward +inf),
    saturating. Mirrors mean in src/lib.fe."""
    acc = 0
    for v in vals:
        acc = clamp_i64(acc + v)
    n = len(vals)
    q = trunc_div(acc, n)
    rem = acc - q * n
    # Round to nearest, ties toward +inf (same convention as sat_mul).
    # Positive: round up iff 2*rem >= n, i.e. rem >= (n+1)//2.
    # Negative: round down iff 2*(-rem) > n, i.e. -rem >= n//2 + 1
    # (a tie at exactly -0.5 rounds toward +inf, i.e. stays at q).
    if acc >= 0:
        if rem >= (n + 1) // 2:
            q += 1
    else:
        if rem <= -(n // 2 + 1):
            q -= 1
    return clamp_i32(q)


def rms_py(vals):
    """RMS of Q15 raws: sqrt of the mean of squares, rounded to nearest.
    Mirrors rms in src/lib.fe."""
    acc = 0
    for v in vals:
        acc = clamp_i64(acc + v * v)
    m = acc // len(vals)  # acc >= 0, trunc == floor
    return clamp_i32(isqrt_round(m))


def biquad_coeffs_lowpass_01():
    """2nd-order Butterworth lowpass, fc = 0.1 * fs, Q15 raw coeffs
    (b0, b1, b2, a1, a2). Standard RBJ/Butterworth design."""
    omega = 2.0 * math.pi * 0.1
    k = math.tan(omega / 2.0)
    norm = 1.0 / (1.0 + math.sqrt(2.0) * k + k * k)
    b0 = k * k * norm
    b1 = 2.0 * b0
    b2 = b0
    a1 = 2.0 * (k * k - 1.0) * norm
    a2 = (1.0 - math.sqrt(2.0) * k + k * k) * norm
    return [int(round(c * 32768)) for c in (b0, b1, b2, a1, a2)]


def biquad_step(state, coeffs, x):
    """One transposed-direct-form-II step, mirroring Biquad::step.
    state = (z1, z2) Q15 raws; returns ((z1', z2'), y)."""
    b0, b1, b2, a1, a2 = coeffs
    z1, z2 = state
    y = sat_add_i32(sat_mul(b0, x, 15), z1)
    z1n = sat_add_i32(sat_sub_i32(sat_mul(b1, x, 15), sat_mul(a1, y, 15)), z2)
    z2n = sat_sub_i32(sat_mul(b2, x, 15), sat_mul(a2, y, 15))
    return (z1n, z2n), y


def biquad_run(coeffs, xs):
    state = (0, 0)
    ys = []
    for x in xs:
        state, y = biquad_step(state, coeffs, x)
        ys.append(y)
    return ys


def fir_step(taps, delay, x):
    """One FIR step, mirroring Fir::step: shift the delay line, then dot.
    Returns (new delay, y)."""
    d = [x] + delay[:-1]
    return d, dot(taps, d, 15)


def fir_run(taps, xs):
    delay = [0] * len(taps)
    ys = []
    for x in xs:
        delay, y = fir_step(taps, delay, x)
        ys.append(y)
    return ys


def sat_neg(x):
    """Saturating negation of an i32 raw, mirroring Q::sat_neg."""
    return I32_MAX if x == I32_MIN else -x


def twiddles(n):
    """FFT twiddle factors W_n^k = exp(-2*pi*i*k/n), Q15 raw.
    Real part: round(cos(2*pi*k/n) * 32768); imaginary part:
    round(-sin(2*pi*k/n) * 32768). Mirrors the TWIDDLE tables in
    src/lib.fe (Q15::one() is 32768, so cos(0) is exactly 32768)."""
    re = [clamp_i32(int(round(math.cos(2.0 * math.pi * k / n) * 32768)))
          for k in range(n)]
    im = [clamp_i32(int(round(-math.sin(2.0 * math.pi * k / n) * 32768)))
          for k in range(n)]
    return re, im


def bitrev(n):
    """Bit-reversal permutation for a radix-2 FFT of size n."""
    log = n.bit_length() - 1
    out = []
    for i in range(n):
        r = 0
        for b in range(log):
            if (i >> b) & 1:
                r |= 1 << (log - 1 - b)
        out.append(r)
    return out


def cplx_mul(a, b):
    """Complex Q15 multiply, mirroring Cplx::mul: (ac-bd) + (ad+bc)i,
    each product a sat_mul, each sum/difference saturating."""
    re = sat_sub_i32(sat_mul(a[0], b[0], 15), sat_mul(a[1], b[1], 15))
    im = sat_add_i32(sat_mul(a[0], b[1], 15), sat_mul(a[1], b[0], 15))
    return (re, im)


def fft_mirror(xs, tw_re, tw_im, rev):
    """Iterative radix-2 DIT FFT, mirroring fft_core in src/lib.fe
    step for step: bit-reversal permutation, then log2(n) stages of
    butterflies with the same twiddle indexing (j * n/len) and the
    same saturating adds/subs."""
    n = len(xs)
    a = [xs[rev[i]] for i in range(n)]
    length = 2
    while length <= n:
        half = length // 2
        step = n // length
        base = 0
        while base < n:
            for j in range(half):
                w = (tw_re[j * step], tw_im[j * step])
                t = cplx_mul(w, a[base + j + half])
                u = a[base + j]
                a[base + j] = (sat_add_i32(u[0], t[0]),
                               sat_add_i32(u[1], t[1]))
                a[base + j + half] = (sat_sub_i32(u[0], t[0]),
                                      sat_sub_i32(u[1], t[1]))
            base += length
        length *= 2
    return a


def dft_float(xs):
    """Reference DFT in floating point (independent of the butterfly
    algorithm): X[k] = sum_n x[n] * exp(-2*pi*i*k*n/N), rounded to
    nearest integer. The mirrored FFT must agree with this within a
    small tolerance (a few ulp per component, growing with log2 N);
    a wrong algorithm would be off by orders of magnitude more."""
    import cmath
    n = len(xs)
    out = []
    for k in range(n):
        s = 0j
        for t, (re, im) in enumerate(xs):
            s += complex(re, im) * cmath.exp(-2j * math.pi * k * t / n)
        out.append((int(round(s.real)), int(round(s.imag))))
    return out


def div_pow2_round(x, s):
    """Divide an i32 raw by 2**s, rounding to nearest (ties toward
    +inf), saturating. Mirrors div_pow2_round in src/lib.fe (the
    inverse-FFT scale by 1/N)."""
    if s <= 0:
        return clamp_i32(x)
    return clamp_i32((x + (1 << (s - 1))) >> s)


def ifft_mirror(xs, tw_re, tw_im, rev):
    """Inverse FFT by the conjugate trick, mirroring ifft_core:
    conj -> fft -> conj -> scale by 1/N with div_pow2_round."""
    n = len(xs)
    log = n.bit_length() - 1
    conj = [(re, sat_neg(im)) for re, im in xs]
    f = fft_mirror(conj, tw_re, tw_im, rev)
    return [(div_pow2_round(re, log), div_pow2_round(sat_neg(im), log))
            for re, im in f]


def cplx_magnitude(re, im):
    """|re + i*im| in Q15, mirroring Cplx::magnitude: sqrt of the
    saturating sum of the squared components."""
    return q_sqrt(sat_add_i32(sat_mul(re, re, 15), sat_mul(im, im, 15)), 15)


def sum_sq(xs):
    """Sum of |x|^2 over a complex frame, exact Python ints."""
    return sum(re * re + im * im for re, im in xs)


def main():
    f15 = 15
    v = {}

    # --- conversions -----------------------------------------------------
    v["from_int"] = {
        "q15_2": from_int(2, 15),
        "q15_neg3": from_int(-3, 15),
        "q15_70000_clamp": from_int(70000, 15),   # > 65535 -> MAX
        "q15_neg70000_clamp": from_int(-70000, 15),
        "q31_1_clamp": from_int(1, 31),           # 1.0 not representable -> MAX
        "q31_neg2_clamp": from_int(-2, 31),
        "q7_100": from_int(100, 7),
    }
    v["from_ratio"] = {
        "q15_3_4": from_ratio(3, 4, 15),          # 0.75 -> 24576
        "q15_neg1_2": from_ratio(-1, 2, 15),      # -0.5 -> -16384
        "q15_1_3": from_ratio(1, 3, 15),          # 0.333.. -> 10923 (rounded)
        "q15_2_3": from_ratio(2, 3, 15),          # 0.666.. -> 21845 (rounded)
        "q15_1_0": from_ratio(1, 0, 15),          # div by zero -> MAX
        "q15_neg1_0": from_ratio(-1, 0, 15),      # -> MIN
        "q15_0_0": from_ratio(0, 0, 15),          # -> 0
        # Negative denominators (v0.4.1 audit: the old formula rounded
        # toward zero here instead of to nearest, ties away from zero).
        "q15_11_neg4": from_ratio(11, -4, 15),    # -2.75 -> -90112
        "q15_neg11_neg4": from_ratio(-11, -4, 15),  # 2.75 -> 90112
        "q15_10_neg4": from_ratio(10, -4, 15),    # -2.5, tie -> -81920
        "q15_neg20_neg7": from_ratio(-20, -7, 15),  # 2.857.. -> 93623
        "q15_1_neg3": from_ratio(1, -3, 15),      # -0.333.. -> -10923
    }

    # --- saturating add/sub ----------------------------------------------
    v["sat_add"] = {
        "max_plus_1": clamp_i32(I32_MAX + 1),
        "min_plus_neg1": clamp_i32(I32_MIN - 1),
        "neg_plus_pos": -100 + 250,
    }
    v["sat_sub"] = {
        "min_minus_1": clamp_i32(I32_MIN - 1),
        "max_minus_neg1": clamp_i32(I32_MAX + 1),
    }

    # --- saturating mul ----------------------------------------------------
    v["sat_mul"] = {
        "half_times_half": sat_mul(16384, 16384, 15),       # 0.25 -> 8192
        "neg_half_times_half": sat_mul(-16384, 16384, 15),  # -0.25 -> -8192
        "one_times_one": sat_mul(32768, 32768, 15),         # 1.0 * 1.0 = 1.0 -> 32768
        "max_times_max": sat_mul(I32_MAX, I32_MAX, 15),     # clamps to MAX
        "min_times_min": sat_mul(I32_MIN, I32_MIN, 15),     # clamps to MAX (pos)
        "min_times_one": sat_mul(I32_MIN, 32768, 15),       # -65536 -> MIN
        "tie_up": sat_mul(1, 16384, 15),                    # +0.5 ulp tie -> 1
        "tie_neg": sat_mul(-1, 16384, 15),                  # -0.5 ulp tie -> 0
        "tie_neg3": sat_mul(-3, 16384, 15),                 # -1.5 ulp -> -1
        "third": sat_mul(from_ratio(1, 3, 15), 32768, 15),  # (1/3)*1.0
        "q31_near1": sat_mul(I32_MAX, I32_MAX, 31),
    }

    # --- saturating div ----------------------------------------------------
    v["sat_div"] = {
        "one_by_two": sat_div(32768, 65536, 15),   # 1.0/2.0 = 0.5 -> 16384
        "three_by_two": sat_div(98304, 65536, 15),  # 3.0/2.0 = 1.5 -> 49152
        "one_by_zero": sat_div(32768, 0, 15),      # -> MAX
        "neg_one_by_zero": sat_div(-32768, 0, 15),  # -> MIN
        "zero_by_zero": sat_div(0, 0, 15),          # -> 0
        "max_by_half": sat_div(I32_MAX, 16384, 15),  # huge -> clamps MAX
        "neg_one_by_three": sat_div(-32768, 98304, 15),  # -1.0/3.0 trunc -> -10922
    }

    # --- scale_int / shifts -------------------------------------------------
    v["scale_int"] = {
        "half_times_3": clamp_i32(16384 * 3),       # 0.5*3 = 1.5 -> 49152
        "max_times_2": clamp_i32(I32_MAX * 2),     # -> MAX
    }
    v["shr"] = {
        "half_shr1": 16384 >> 1,                   # 0.25 -> 8192
        "neg_shr1": -16384 >> 1,                   # -0.25 -> -8192
        "pos_shr40": 0,                            # n>=31 -> 0
        "neg_shr40": -1,                           # n>=31 -> -1
    }
    v["shl_sat"] = {
        "half_shl1": clamp_i32(16384 << 1),        # 1.0 -> 32768
        "max_shl1": clamp_i32(I32_MAX << 1),       # -> MAX
        "neg_shl40": I32_MIN,                      # n>=32, negative -> MIN
        "zero_shl40": 0,
    }

    # --- abs / neg ----------------------------------------------------------
    v["abs_neg"] = {
        "abs_min": I32_MAX,                        # saturating abs(MIN)
        "abs_neg5": 5,
        "neg_min": I32_MAX,                        # saturating neg(MIN)
        "neg_5": -5,
    }

    # --- Hann tables ---------------------------------------------------------
    tables = {n: hann(n) for n in (16, 32, 64)}
    v["hann"] = {
        f"hann{n}": t for n, t in tables.items()
    }
    v["hann_props"] = {
        "hann16_sum": sum(tables[16]),
        "hann32_sum": sum(tables[32]),
        "hann64_sum": sum(tables[64]),
        "hann16_peak": max(tables[16]),
        "hann32_peak": max(tables[32]),
        "hann64_peak": max(tables[64]),
    }

    # --- dot products --------------------------------------------------------
    rng = random.Random(42)
    a4 = [rng.randint(-30000, 30000) for _ in range(4)]
    b4 = [rng.randint(-30000, 30000) for _ in range(4)]
    a16 = [rng.randint(-30000, 30000) for _ in range(16)]
    b16 = [rng.randint(-30000, 30000) for _ in range(16)]
    v["dot"] = {
        "a4": a4, "b4": b4, "dot4": dot(a4, b4, 15),
        "a16": a16, "b16": b16, "dot16": dot(a16, b16, 15),
        "unit4": dot([16384] * 4, [16384] * 4, 15),  # 4 * 0.25 = 1.0 -> 32768
    }

    # --- demo: sine burst -> Hann window -> 0.5 gain -> energy ---------------
    burst = [
        int(round(math.sin(2.0 * math.pi * 2 * n / 16.0) * 32767 * 0.9))
        for n in range(16)
    ]
    h16 = tables[16]
    windowed = [sat_mul(x, w, 15) for x, w in zip(burst, h16)]
    gained = [sat_mul(x, 16384, 15) for x in windowed]  # x 0.5
    energy = dot(gained, gained, 15)
    v["demo"] = {
        "burst": burst,
        "windowed": windowed,
        "gained": gained,
        "energy": energy,
        "windowed_sum": sum(windowed),
    }

    # --- dot saturation (v0.2.0 regression) --------------------------------
    # 3 max-magnitude elements used to revert (i64 accumulator overflow);
    # saturating accumulation must return MAX / MIN instead.
    v["dot_sat"] = {
        "max3": dot([I32_MAX] * 3, [I32_MAX] * 3, 15),
        "min3": dot([I32_MIN] * 3, [I32_MIN] * 3, 15),
        "max_min3": dot([I32_MAX] * 3, [I32_MIN] * 3, 15),
        "min_max3": dot([I32_MIN] * 3, [I32_MAX] * 3, 15),
        "max8": dot([I32_MAX] * 8, [I32_MAX] * 8, 15),
        "min8": dot([I32_MIN] * 8, [I32_MIN] * 8, 15),
        "mixed": dot([I32_MAX, I32_MIN, 16384, -16384],
                     [16384, -16384, I32_MAX, I32_MIN], 15),
    }

    # --- format conversion ---------------------------------------------------
    v["convert_q"] = {
        "q15_half_to_q31": convert_q(16384, 15, 31),
        "q15_half_to_q7": convert_q(16384, 15, 7),
        "q7_one_to_q15": convert_q(128, 7, 15),
        "q31_max_to_q15": convert_q(I32_MAX, 31, 15),
        "q31_min_to_q15": convert_q(I32_MIN, 31, 15),
        "q15_max_to_q7": convert_q(I32_MAX, 15, 7),
        "q15_one_to_q7": convert_q(32768, 15, 7),
        "q15_raw1_to_q7": convert_q(1, 15, 7),
        "q15_raw128_to_q7": convert_q(128, 15, 7),
        "q15_neg128_to_q7": convert_q(-128, 15, 7),
        "q15_neg129_to_q7": convert_q(-129, 15, 7),
        "q15_neg1_to_q7": convert_q(-1, 15, 7),
        "q0_5_to_q15": convert_q(5, 0, 15),
        "q15_5_to_q0": convert_q(163840, 15, 0),
        "q15_neg_to_q0": convert_q(-49152, 15, 0),
        "q0_1_to_q31": convert_q(1, 0, 31),
        "q31_oneish_to_q0": convert_q(I32_MAX, 31, 0),
        "q31_neg_to_q0": convert_q(I32_MIN, 31, 0),
        "q15_to_q15": convert_q(-12345, 15, 15),
    }

    # --- Hamming / Blackman tables --------------------------------------------
    ham = {n: hamming(n) for n in (16, 32, 64)}
    blk = {n: blackman(n) for n in (16, 32, 64)}
    v["hamming"] = {f"hamming{n}": t for n, t in ham.items()}
    v["blackman"] = {f"blackman{n}": t for n, t in blk.items()}
    v["window_props"] = {
        "hamming16_sum": sum(ham[16]), "hamming32_sum": sum(ham[32]),
        "hamming64_sum": sum(ham[64]),
        "hamming16_peak": max(ham[16]), "hamming32_peak": max(ham[32]),
        "hamming64_peak": max(ham[64]),
        "hamming16_end": ham[16][0], "hamming16_mid": ham[16][8],
        "blackman16_sum": sum(blk[16]), "blackman32_sum": sum(blk[32]),
        "blackman64_sum": sum(blk[64]),
        "blackman16_peak": max(blk[16]), "blackman32_peak": max(blk[32]),
        "blackman64_peak": max(blk[64]),
        "blackman16_end": blk[16][0], "blackman16_mid": blk[16][8],
    }

    # --- mean / rms / energy ---------------------------------------------------
    demo_frame = burst  # the demo sine burst, reused as a frame
    const_frame = [16384] * 16
    v["frame_stats"] = {
        "burst_mean": mean_py(demo_frame),
        "burst_rms": rms_py(demo_frame),
        "burst_energy": dot(demo_frame, demo_frame, 15),
        "const_mean": mean_py(const_frame),
        "const_rms": rms_py(const_frame),
        "const_energy": dot(const_frame, const_frame, 15),
        "mean_1_2": mean_py([1, 2]),
        "mean_neg1_neg2": mean_py([-1, -2]),
        "mean_max_pair": mean_py([I32_MAX, I32_MAX]),
        "mean_min_pair": mean_py([I32_MIN, I32_MIN]),
        "mean_max_min": mean_py([I32_MAX, I32_MIN]),
        "rms_half_zero": rms_py([16384, 0, 0, 0]),
        "rms_full": rms_py([32768, 32768, 32768, 32768]),
        "rms_min": rms_py([I32_MIN]),
        "rms_zero": rms_py([0, 0, 0, 0]),
        "energy_unit": dot([16384] * 4, [16384] * 4, 15),
    }

    # --- sqrt / reciprocal ------------------------------------------------------
    v["sqrt"] = {
        "q15_one": q_sqrt(32768, 15),
        "q15_quarter": q_sqrt(8192, 15),
        "q15_half": q_sqrt(16384, 15),
        "q15_two": q_sqrt(65536, 15),
        "q15_four": q_sqrt(131072, 15),
        "q15_zero": q_sqrt(0, 15),
        "q15_neg": q_sqrt(-32768, 15),
        "q15_009": q_sqrt(2949, 15),
        "q15_max": q_sqrt(I32_MAX, 15),
        "q31_max": q_sqrt(I32_MAX, 31),
        "q31_half": q_sqrt(1073741824, 31),
        "q7_four": q_sqrt(512, 7),
        "q7_two": q_sqrt(256, 7),
        "q0_16": q_sqrt(16, 0),
        "q0_2": q_sqrt(2, 0),
        "q0_0": q_sqrt(0, 0),
    }
    v["reciprocal"] = {
        "q15_one": q_reciprocal(32768, 15),
        "q15_half": q_reciprocal(16384, 15),
        "q15_two": q_reciprocal(65536, 15),
        "q15_zero": q_reciprocal(0, 15),
        "q15_neg2": q_reciprocal(-65536, 15),
        "q15_third": q_reciprocal(10923, 15),
        "q15_three": q_reciprocal(98304, 15),
        "q15_quarter": q_reciprocal(8192, 15),
        "q15_raw1": q_reciprocal(1, 15),
        "q15_raw2": q_reciprocal(2, 15),
        "q15_raw7": q_reciprocal(7, 15),
        "q15_max": q_reciprocal(I32_MAX, 15),
        "q15_min": q_reciprocal(I32_MIN, 15),
        "q15_neg1": q_reciprocal(-1, 15),
        "q7_two": q_reciprocal(256, 7),
        "q7_half": q_reciprocal(64, 7),
        "q7_one": q_reciprocal(128, 7),
        "q31_max": q_reciprocal(I32_MAX, 31),
        "q31_min": q_reciprocal(I32_MIN, 31),
        "q31_half": q_reciprocal(1073741824, 31),
        "q0_5": q_reciprocal(5, 0),
        "q0_1": q_reciprocal(1, 0),
        "q0_neg1": q_reciprocal(-1, 0),
        "q0_0": q_reciprocal(0, 0),
    }

    # --- biquad (Butterworth lowpass, fc = 0.1 fs) ------------------------------
    lp = biquad_coeffs_lowpass_01()
    v["biquad"] = {
        "lp_coeffs": lp,
        "impulse8": biquad_run(lp, [32768] + [0] * 7),
        "step16": biquad_run(lp, [32768] * 16),
        "step16_last": biquad_run(lp, [32768] * 16)[-1],
        "gain_half_impulse": biquad_run([16384, 0, 0, 0, 0], [16384, 0, 0, 0]),
        "gain_half_step": biquad_run([16384, 0, 0, 0, 0], [16384] * 4),
    }

    # --- FIR ---------------------------------------------------------------------
    rng7 = random.Random(7)
    fir_taps = [rng7.randint(-30000, 30000) for _ in range(8)]
    ma4 = [8192, 8192, 8192, 8192]
    v["fir"] = {
        "taps8": fir_taps,
        "taps8_impulse": fir_run(fir_taps, [32768] + [0] * 9),
        "taps8_step": fir_run(fir_taps, [32768] * 10),
        "ma4_impulse": fir_run(ma4, [32768] + [0] * 5),
        "ma4_step": fir_run(ma4, [32768] * 6),
        "ma4_const": fir_run(ma4, [16384] * 6),
    }

    # --- FFT (v0.3.0) ---------------------------------------------------------
    # Twiddle + bit-reversal tables, and FFT vectors for 16/32/64.
    # Cross-checks (fail loudly at generation time if violated):
    #   * the mirrored butterfly agrees with an independent float DFT
    #     within 16 ulp per component (observed max: 11 at N=64);
    #   * forward-then-inverse round-trips return the input within 2 ulp
    #     (observed max: 1);
    #   * Parseval: sum |X|^2 == N * sum |x|^2 within the bound below
    #     (observed diffs printed; bound set with headroom).
    fft_tw = {n: twiddles(n) for n in (16, 32, 64)}
    fft_rev = {n: bitrev(n) for n in (16, 32, 64)}
    v["fft_twiddle"] = {
        f"tw{n}_re": fft_tw[n][0] for n in (16, 32, 64)
    }
    v["fft_twiddle"].update({
        f"tw{n}_im": fft_tw[n][1] for n in (16, 32, 64)
    })
    v["fft_bitrev"] = {f"bitrev{n}": fft_rev[n] for n in (16, 32, 64)}

    fft_inputs = {}
    for n in (16, 32, 64):
        impulse = [(32768, 0)] + [(0, 0)] * (n - 1)
        dc = [(32768, 0)] * n
        sine1 = [(int(round(math.sin(2.0 * math.pi * t / n) * 32767 * 0.9)), 0)
                 for t in range(n)]
        sine2 = [(int(round(math.sin(2.0 * math.pi * 2 * t / n) * 32767 * 0.9)), 0)
                 for t in range(n)]
        rng99 = random.Random(99)
        rnd = [(rng99.randint(-30000, 30000), rng99.randint(-30000, 30000))
               for _ in range(n)]
        fft_inputs[n] = {"impulse": impulse, "dc": dc, "sine1": sine1,
                         "sine2": sine2, "random": rnd}

    v["fft"] = {}
    parseval_diffs = {}
    for n in (16, 32, 64):
        tw_re, tw_im = fft_tw[n]
        rev = fft_rev[n]
        for name, xs in fft_inputs[n].items():
            m = fft_mirror(xs, tw_re, tw_im, rev)
            d = dft_float(xs)
            err = max(max(abs(mr - dr), abs(mi - di))
                      for (mr, mi), (dr, di) in zip(m, d))
            assert err <= 16, f"fft n={n} {name}: mirror vs float DFT err {err} > 16"
            rt = ifft_mirror(m, tw_re, tw_im, rev)
            rt_err = max(max(abs(a - c), abs(b - dd))
                         for (a, b), (c, dd) in zip(rt, xs))
            assert rt_err <= 2, f"fft n={n} {name}: round-trip err {rt_err} > 2"
            tsum = sum_sq(xs)
            fsum = sum_sq(m)
            pdiff = abs(fsum - n * tsum)
            parseval_diffs[f"n{n}_{name}"] = pdiff
            v["fft"][f"n{n}_{name}"] = {
                "input_re": [re for re, _ in xs],
                "input_im": [im for _, im in xs],
                "out_re": [re for re, _ in m],
                "out_im": [im for _, im in m],
                "magnitude": [cplx_magnitude(re, im) for re, im in m],
                "roundtrip_re": [re for re, _ in rt],
                "roundtrip_im": [im for _, im in rt],
                "roundtrip_err": rt_err,
                "dft_err": err,
                "time_sq_sum": tsum,
                "freq_sq_sum": fsum,
                "parseval_diff": pdiff,
            }
    print("fft mirror-vs-DFT max err:",
          max(v["fft"][k]["dft_err"] for k in v["fft"]))
    print("fft round-trip max err:",
          max(v["fft"][k]["roundtrip_err"] for k in v["fft"]))
    print("fft parseval diffs:", parseval_diffs)

    # --- Taylor-Green slice (Rank 1, Clay mapping) -------------------
    # 1D slice of the Taylor-Green vorticity field w = 2 sin x sin y at
    # sin y = 1: x[n] = 2 sin(2*pi*n/N), amplitude exactly 2.0 in Q15.
    # Its DFT has exactly two non-zero bins: bin 1 = -i*N, bin N-1 =
    # +i*N (Q15 values), i.e. im raws -/+ N*32768. Every other bin is
    # zero in exact arithmetic; the mirrored FFT reports the small
    # residuals from input/twiddle quantization. Cross-check bounds
    # are looser than the fixed FFT vectors (16 / 2 ulp) because the
    # error scales with input amplitude (2.0 here vs <= 1.0 there) and
    # with N: observed mirror-vs-DFT 3/9/27 and round-trip 2/2/3 at
    # N = 16/32/64, i.e. ~1e-5 relative on the peak bin at N=64.
    # Bounds below (32 / 4 ulp) carry headroom over those observations.
    v["taylor_green"] = {}
    for n in (16, 32, 64):
        xs = [(int(round(2.0 * math.sin(2.0 * math.pi * t / n) * 32768)), 0)
              for t in range(n)]
        tw_re, tw_im = fft_tw[n]
        rev = fft_rev[n]
        m = fft_mirror(xs, tw_re, tw_im, rev)
        d = dft_float(xs)
        err = max(max(abs(mr - dr), abs(mi - di))
                  for (mr, mi), (dr, di) in zip(m, d))
        assert err <= 32, f"taylor-green n={n}: mirror vs float DFT err {err} > 32"
        rt = ifft_mirror(m, tw_re, tw_im, rev)
        rt_err = max(max(abs(a - c), abs(b - dd))
                     for (a, b), (c, dd) in zip(rt, xs))
        assert rt_err <= 4, f"taylor-green n={n}: round-trip err {rt_err} > 4"
        peak = n * 32768
        assert abs(m[1][0]) <= 32 and abs(m[1][1] + peak) <= 32, \
            f"taylor-green n={n}: bin 1 = {m[1]}, expected (0, {-peak})"
        assert abs(m[n - 1][0]) <= 32 and abs(m[n - 1][1] - peak) <= 32, \
            f"taylor-green n={n}: bin {n - 1} = {m[n - 1]}, expected (0, {peak})"
        other_max = max(max(abs(m[k][0]), abs(m[k][1]))
                        for k in range(n) if k not in (1, n - 1))
        assert other_max <= 32, f"taylor-green n={n}: non-peak bin {other_max} > 32"
        v["taylor_green"][f"n{n}"] = {
            "input": [re for re, _ in xs],
            "out_re": [re for re, _ in m],
            "out_im": [im for _, im in m],
            "magnitude": [cplx_magnitude(re, im) for re, im in m],
            "roundtrip": [re for re, _ in rt],
            "roundtrip_im": [im for _, im in rt],
            "roundtrip_err": rt_err,
            "dft_err": err,
            "peak_raw": peak,
            "other_max": other_max,
        }
    print("taylor-green:", {k: (v["taylor_green"][k]["dft_err"],
                                v["taylor_green"][k]["roundtrip_err"],
                                v["taylor_green"][k]["other_max"])
                           for k in v["taylor_green"]})

    # --- Cplx unit values (v0.3.0) ---------------------------------------------
    half_c = (16384, 16384)  # 0.5 + 0.5i
    v["cplx"] = {
        "half_sq": cplx_mul(half_c, half_c),          # (0.5+0.5i)^2 = 0 + 0.5i
        "half_conj": (16384, sat_neg(16384)),
        "mag_0608": cplx_magnitude(19661, 26214),     # |0.6 + 0.8i| = 1.0
        "mag_sq_half": sat_add_i32(sat_mul(16384, 16384, 15),
                                   sat_mul(16384, 16384, 15)),  # 0.5
        "mag_one": cplx_magnitude(32768, 0),
        "mag_zero": cplx_magnitude(0, 0),
        "add": (sat_add_i32(16384, 8192), sat_add_i32(-16384, 16384)),
        "sub": (sat_sub_i32(16384, 8192), sat_sub_i32(-16384, 16384)),
        "neg_conj_min": sat_neg(I32_MIN),             # conj of MIN im -> MAX
    }

    # --- demo signal chain (v0.3.0) -------------------------------------------
    # Sine burst -> Hann window -> 4-tap moving-average FIR ->
    # Butterworth lowpass biquad -> RMS/energy -> FFT magnitude.
    ma4_chain = [8192, 8192, 8192, 8192]
    fir_out = fir_run(ma4_chain, windowed)
    bq_out = biquad_run(lp, fir_out)
    spec = fft_mirror([(x, 0) for x in bq_out],
                      fft_tw[16][0], fft_tw[16][1], fft_rev[16])
    spec_mag = [cplx_magnitude(re, im) for re, im in spec]
    v["demo_chain"] = {
        "fir_out": fir_out,
        "bq_out": bq_out,
        "rms": rms_py(bq_out),
        "energy": dot(bq_out, bq_out, 15),
        "fft_re": [re for re, _ in spec],
        "fft_im": [im for _, im in spec],
        "fft_mag": spec_mag,
    }

    # --- randomized differential tests (v0.4.0) ---------------------------
    # Seeded deterministic RNG: seed 20261120, consumed in the fixed
    # order dot -> fir -> biquad -> fft -> scalars -> demo-contract
    # frame, so re-running this script regenerates byte-identical
    # vectors. Differential coverage on top of the fixed vectors above:
    # random inputs through the same functions, with every expected
    # value computed here, never copied from Fe output.
    drng = random.Random(20261120)

    # dot: 4 random 8-length pairs + 2 random 16-length pairs.
    dot_rand = []
    for _ in range(4):
        a = [drng.randint(-30000, 30000) for _ in range(8)]
        b = [drng.randint(-30000, 30000) for _ in range(8)]
        dot_rand.append({"a": a, "b": b, "out": dot(a, b, 15)})
    for _ in range(2):
        a = [drng.randint(-30000, 30000) for _ in range(16)]
        b = [drng.randint(-30000, 30000) for _ in range(16)]
        dot_rand.append({"a": a, "b": b, "out": dot(a, b, 15)})

    # FIR: random 8 taps + random 16-sample input.
    fir_taps_r = [drng.randint(-30000, 30000) for _ in range(8)]
    fir_in_r = [drng.randint(-30000, 30000) for _ in range(16)]
    fir_out_r = fir_run(fir_taps_r, fir_in_r)

    # biquad: random 16-sample input through the documented lowpass.
    bq_in_r = [drng.randint(-30000, 30000) for _ in range(16)]
    bq_out_r = biquad_run(lp, bq_in_r)

    # FFT: random complex frames, 16 and 32 points. Same cross-checks
    # as the fixed FFT vectors: mirrored butterfly vs independent
    # float DFT (<= 16 ulp), forward-then-inverse round-trip (<= 2 ulp).
    fft_rand = {}
    for n in (16, 32):
        xs = [(drng.randint(-30000, 30000), drng.randint(-30000, 30000))
              for _ in range(n)]
        tw_re, tw_im = fft_tw[n]
        rev = fft_rev[n]
        m = fft_mirror(xs, tw_re, tw_im, rev)
        d = dft_float(xs)
        err = max(max(abs(mr - dr), abs(mi - di))
                  for (mr, mi), (dr, di) in zip(m, d))
        assert err <= 16, f"fft random n={n}: mirror vs float DFT err {err} > 16"
        rt = ifft_mirror(m, tw_re, tw_im, rev)
        rt_err = max(max(abs(a - c), abs(b - dd))
                     for (a, b), (c, dd) in zip(rt, xs))
        assert rt_err <= 2, f"fft random n={n}: round-trip err {rt_err} > 2"
        fft_rand[f"n{n}"] = {
            "in_re": [re for re, _ in xs],
            "in_im": [im for _, im in xs],
            "out_re": [re for re, _ in m],
            "out_im": [im for _, im in m],
            "magnitude": [cplx_magnitude(re, im) for re, im in m],
            "roundtrip_re": [re for re, _ in rt],
            "roundtrip_im": [im for _, im in rt],
            "dft_err": err,
            "roundtrip_err": rt_err,
            "time_sq_sum": sum_sq(xs),
            "freq_sq_sum": sum_sq(m),
        }

    # scalars: 16 random raws for sqrt/reciprocal (10 in audio range,
    # 3 near zero, 2 near the i32 extremes, 1 exactly), plus convert_q
    # over 6 format pairs for the first 10 raws.
    scal_r = []
    for _ in range(10):
        scal_r.append(drng.randint(-30000, 30000))
    for _ in range(3):
        scal_r.append(drng.randint(-100, 100))
    scal_r.append(drng.randint(I32_MIN // 2, I32_MAX // 2))
    scal_r.append(drng.randint(-I32_MAX, I32_MAX))
    scal_r.append(1)
    scal_pairs = [(15, 7), (15, 31), (7, 15), (31, 15), (15, 0), (0, 15)]
    conv_r = []
    for raw in scal_r[:10]:
        for frm, to in scal_pairs:
            conv_r.append({"raw": raw, "from": frm, "to": to, "out": convert_q(raw, frm, to)})

    # demo-contract frame: a random 16-sample frame run through the full
    # demo chain (Hann window -> 4-tap MA FIR -> Butterworth biquad ->
    # RMS/energy -> FFT magnitude), for the DspDemo contract test.
    dc_frame = [drng.randint(-30000, 30000) for _ in range(16)]
    dc_windowed = [sat_mul(x, w, 15) for x, w in zip(dc_frame, tables[16])]
    dc_fir = fir_run(ma4_chain, dc_windowed)
    dc_bq = biquad_run(lp, dc_fir)
    dc_spec = fft_mirror([(x, 0) for x in dc_bq],
                         fft_tw[16][0], fft_tw[16][1], fft_rev[16])
    dc_mag = [cplx_magnitude(re, im) for re, im in dc_spec]
    dc_peak = max(dc_mag)
    dc_bin = dc_mag.index(dc_peak)

    v["differential"] = {
        "seed": 20261120,
        "dot": dot_rand,
        "fir": {"taps": fir_taps_r, "input": fir_in_r, "output": fir_out_r},
        "biquad": {"input": bq_in_r, "output": bq_out_r},
        "fft": fft_rand,
        "scalars": scal_r,
        "sqrt": [q_sqrt(r, 15) for r in scal_r],
        "reciprocal": [q_reciprocal(r, 15) for r in scal_r],
        "convert": conv_r,
    }
    v["demo_contract"] = {
        "frame": dc_frame,
        "rms": rms_py(dc_bq),
        "energy": dot(dc_bq, dc_bq, 15),
        "fft_mag": dc_mag,
        "fft_peak": dc_peak,
        "fft_peak_bin": dc_bin,
    }

    with open("tools/vectors.json", "w") as fh:
        json.dump(v, fh, indent=2)

    with open("tools/hann_tables.txt", "w") as fh:
        for n in (16, 32, 64):
            t = tables[n]
            fh.write(f"HANN{n} ({len(t)} values):\n")
            fh.write(", ".join(str(x) for x in t) + "\n\n")

    with open("tools/window_tables.txt", "w") as fh:
        for name, tbl in (("HAMMING", ham), ("BLACKMAN", blk)):
            for n in (16, 32, 64):
                t = tbl[n]
                fh.write(f"{name}{n} ({len(t)} values):\n")
                fh.write(", ".join(str(x) for x in t) + "\n\n")

    with open("tools/fft_tables.txt", "w") as fh:
        for n in (16, 32, 64):
            tw_re, tw_im = fft_tw[n]
            fh.write(f"TWIDDLE{n}_RE ({len(tw_re)} values):\n")
            fh.write(", ".join(str(x) for x in tw_re) + "\n\n")
            fh.write(f"TWIDDLE{n}_IM ({len(tw_im)} values):\n")
            fh.write(", ".join(str(x) for x in tw_im) + "\n\n")
            fh.write(f"BITREV{n} ({len(fft_rev[n])} values):\n")
            fh.write(", ".join(str(x) for x in fft_rev[n]) + "\n\n")

    print("vectors.json + hann_tables.txt + window_tables.txt + fft_tables.txt written")
    print("demo chain rms:", v["demo_chain"]["rms"],
          "energy:", v["demo_chain"]["energy"])
    print("demo chain fft mag:", v["demo_chain"]["fft_mag"])
    print("differential seed:", v["differential"]["seed"],
          "fft random errs:",
          {k: (v["differential"]["fft"][k]["dft_err"],
               v["differential"]["fft"][k]["roundtrip_err"])
           for k in v["differential"]["fft"]})
    print("demo contract frame rms:", v["demo_contract"]["rms"],
          "energy:", v["demo_contract"]["energy"],
          "fft peak:", v["demo_contract"]["fft_peak"],
          "bin:", v["demo_contract"]["fft_peak_bin"])
    print("hann16 sum:", v["hann_props"]["hann16_sum"],
          "peak:", v["hann_props"]["hann16_peak"])
    print("hann32 sum:", v["hann_props"]["hann32_sum"],
          "peak:", v["hann_props"]["hann32_peak"])
    print("hann64 sum:", v["hann_props"]["hann64_sum"],
          "peak:", v["hann_props"]["hann64_peak"])
    print("dot4:", v["dot"]["dot4"], "dot16:", v["dot"]["dot16"])
    print("demo energy:", energy, "windowed sum:", sum(windowed))
    print("dot_sat max3:", v["dot_sat"]["max3"], "min3:", v["dot_sat"]["min3"],
          "max_min3:", v["dot_sat"]["max_min3"])
    print("biquad lp coeffs:", v["biquad"]["lp_coeffs"])
    print("biquad impulse8:", v["biquad"]["impulse8"])
    print("biquad step16 last:", v["biquad"]["step16_last"])
    print("fir taps8:", v["fir"]["taps8"])
    print("fir taps8 impulse:", v["fir"]["taps8_impulse"])
    print("hamming16 sum:", v["window_props"]["hamming16_sum"],
          "peak:", v["window_props"]["hamming16_peak"],
          "end:", v["window_props"]["hamming16_end"])
    print("blackman16 sum:", v["window_props"]["blackman16_sum"],
          "peak:", v["window_props"]["blackman16_peak"],
          "end:", v["window_props"]["blackman16_end"])


if __name__ == "__main__":
    main()
