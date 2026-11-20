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
  - Hann window: round(0.5 * (1 - cos(2*pi*n/(N-1))) * 32767), Q15 raw.

Outputs:
  - tools/vectors.json : every vector used by the Fe tests.
  - tools/hann_tables.txt : Fe literals for HANN16/32/64.
Run: python3 tools/gen_vectors.py
"""
import json
import math
import random

I32_MAX = 2**31 - 1
I32_MIN = -(2**31)


def clamp_i32(x):
    return max(I32_MIN, min(I32_MAX, x))


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
    # round to nearest, ties away from zero
    half = abs(den) // 2
    if v >= 0:
        return clamp_i32(trunc_div(v + half, den)) if den > 0 else clamp_i32(trunc_div(v - half, den))
    else:
        return clamp_i32(trunc_div(v - half, den)) if den > 0 else clamp_i32(trunc_div(v + half, den))


def hann(n_points):
    return [
        int(round(0.5 * (1.0 - math.cos(2.0 * math.pi * n / (n_points - 1))) * 32767))
        for n in range(n_points)
    ]


def dot(a, b, frac):
    acc = sum(x * y for x, y in zip(a, b))
    if frac >= 1:
        acc = (acc + (1 << (frac - 1))) >> frac
    return clamp_i32(acc)


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

    with open("tools/vectors.json", "w") as fh:
        json.dump(v, fh, indent=2)

    with open("tools/hann_tables.txt", "w") as fh:
        for n in (16, 32, 64):
            t = tables[n]
            fh.write(f"HANN{n} ({len(t)} values):\n")
            fh.write(", ".join(str(x) for x in t) + "\n\n")

    print("vectors.json + hann_tables.txt written")
    print("hann16 sum:", v["hann_props"]["hann16_sum"],
          "peak:", v["hann_props"]["hann16_peak"])
    print("hann32 sum:", v["hann_props"]["hann32_sum"],
          "peak:", v["hann_props"]["hann32_peak"])
    print("hann64 sum:", v["hann_props"]["hann64_sum"],
          "peak:", v["hann_props"]["hann64_peak"])
    print("dot4:", v["dot"]["dot4"], "dot16:", v["dot"]["dot16"])
    print("demo energy:", energy, "windowed sum:", sum(windowed))


if __name__ == "__main__":
    main()
