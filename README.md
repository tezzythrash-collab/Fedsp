# FeDSP — Fixed-Point DSP Math Library for Fe

**One-sentence pitch:** Saturating binary Q-format fixed-point arithmetic for
audio-rate DSP in pure Fe — Q7/Q15/Q31, Hann windows, and dot products, with
every result verified against an independent Python oracle.

**Fe-specific highlight:** the library is built on Fe 26.5 features that do
not exist in Solidity: a generic `Q<const FRAC: i32>` type (one
implementation serves every Q-format), a manual `Copy` impl that makes
fixed-size `Q` arrays usable with repetition and indexing, `i64` widening
with `downcast_saturate()` for overflow-free saturating arithmetic, and pure
functions that run unchanged on both the sonatina and native test backends.

## Problem

Fe's standard library ships **decimal** fixed-point (`core::fixed::Fixed<D>`
= `raw / 10^D`, `Wad` / `Ray`) — built for money and DeFi. Audio and signal
processing need **binary** Q-format (`raw / 2^FRAC`): powers of two map to
shifts, Q15 is the native sample format of audio hardware and codecs, and
DSP code needs saturation (never wrap, never revert) rather than the
checked-arithmetic revert behavior of plain integer math. Before FeDSP,
writing a saturating Q15 gain stage, a Hann window, or a dot product in Fe
meant hand-rolling the widening, rounding, and clamping in every contract.

## Approach

- **One generic type.** `Q<const FRAC: i32>` stores a signed `i32` raw
  value; `Q7`, `Q15`, `Q31` are type aliases. `Q<0>` is plain integers.
- **Saturation everywhere.** Add, sub, mul, div, conversions, shifts, abs,
  and neg all saturate to `i32::MAX` / `i32::MIN` instead of wrapping or
  reverting. Division by zero saturates by sign (or `0` for `0 / 0`).
- **Exact rounding, documented.** Multiplication widens to `i64`, adds half
  an ulp, shifts, and clamps — round-to-nearest, ties toward +inf. Division
  truncates toward zero. `to_int` is an arithmetic shift (floors negatives:
  `-1.5 -> -2`). Every rule is stated in the doc comments and pinned by a
  test.
- **Verified, not asserted.** Every expected value in `src/tests.fe` was
  computed by `tools/gen_vectors.py`, an independent Python oracle using
  exact arbitrary-precision integers — not by copying Fe output. If Fe and
  the oracle ever disagree, the test fails and the library is wrong.

## Architecture

```
fedsp/
  fe.toml                  ingot: fedsp 0.1.0 (Fe 26.5.0+)
  src/lib.fe               the whole library (single file, ~420 lines)
    Q<const FRAC: i32>     conversions, saturating arithmetic, scaling,
                           shifts, clamp/min/max, abs/neg
    Q7 / Q15 / Q31         type aliases
    HANN16 / 32 / 64       Hann window tables, Q15 raw (Python-generated)
    dot<const FRAC, N>     dot product, i64 accumulation, rounded + clamped
    window<const N>        elementwise Hann (or any) windowing of Q15 frames
  src/tests.fe             13 tests, 105 assertions, all oracle-verified
  tools/gen_vectors.py     the Python oracle (vectors + Hann tables)
  tools/vectors.json       oracle output (all expected values)
  .github/workflows/ci.yml  fmt check, fe check, fe test (sonatina + native)
```

Key design facts discovered while building (all verified by compiler probe):

- `core::fixed` is decimal, not Q-format — FeDSP complements it.
- `i32` has `saturating_add` / `saturating_sub` / `saturating_mul` but no
  `checked_add` and no `i32::MAX` / `i32::MIN` constants (bounds are
  literals `2147483647` / `-2147483648`).
- Narrowing casts need `.downcast_saturate()` / `.downcast_truncate()`;
  `as` is only for provably lossless widening.
- `Q` needs a manual `impl Copy for Q<FRAC>` before arrays of `Q` support
  repetition (`[Q15::zero(); N]`) and indexing.
- `Q<31>::one()` cannot be `1.0` (not representable); `one()` saturates to
  `maximum()` instead of silently wrapping to `i32::MIN` (`1 << 31` wraps —
  verified by probe).
- Function arguments are labeled; no turbofish (`Q<15>::from_int(...)`).

## Demo (reproducible CLI example)

The `demo_gain_stage` test is the demo: a 16-sample sine burst (2 cycles,
amplitude 0.9) goes through a Hann window, a 0.5 gain stage, and an energy
measurement — the whole chain in pure Fe, with oracle-verified values at
each step (windowed sum `-2183`, energy `18684`).

```console
$ fe test
PASS  [0.0180s] demo_gain_stage
PASS  [0.0000s] q0_plain_integers
PASS  [0.0000s] sat_add_sub
PASS  [0.0000s] sat_div
PASS  [0.0000s] sat_mul
PASS  [0.0000s] scale_and_shifts
PASS  [0.0226s] constants
PASS  [0.0213s] conversions_from_i64
PASS  [0.0209s] conversions_from_int
PASS  [0.0243s] abs_neg_clamp_min_max
PASS  [0.0181s] dot_products
PASS  [0.0151s] hann_tables
PASS  [0.0231s] conversions_from_ratio

test result: ok. 13 passed; 0 failed
```

The same 13 tests also pass on the native backend (`fe test --backend native`).

## Build & test

Requires Fe **26.5.0** or later (verified with `fe 26.5.0 (3dab95f)`).

```console
# Install the toolchain (Linux x86_64; also works on macOS/Windows via feup)
curl -fsSL https://raw.githubusercontent.com/argotorg/fe/master/feup/feup.sh -o /tmp/feup.sh
bash /tmp/feup.sh --version v26.5.0 --no-modify-path
export PATH="$HOME/.fe/bin:$PATH"

# In the repo root:
fe fmt --check   # formatting
fe check         # type check
fe test          # 13 tests, sonatina backend
fe test --backend native   # 13 tests, native backend
```

`fe build` is intentionally not part of the workflow: FeDSP is a pure
library with no contracts, and `fe build` exits 1 with "No contracts found
to build". `fe check` + `fe test` are the complete build surface.

To regenerate the test vectors and Hann tables:

```console
python3 tools/gen_vectors.py   # rewrites tools/vectors.json + tools/hann_tables.txt
```

## Limitations

- **i32 raw only.** No `i64`-raw or `u32` Q-formats yet; `Q31` tops out
  just below `1.0` (see `one()` above).
- **FRAC range `0..=31`, enforced at compile time.** `Q` carries a
  `where FRAC >= 0, FRAC <= 31` bound (propagated to `dot`); instantiating
  `Q` outside the range is a compile error, not undefined behavior.
- **No FFT, filters, or resampling (yet).** The library covers the
  arithmetic core, windowing, and dot products — the primitives a filter or
  FFT would be built on. Those are natural follow-ups, not included here.
- **No `f32` interop.** Conversion to/from float is out of scope for this
  version; all boundaries are integer.
- **Rounding is fixed.** Mul rounds to nearest (ties toward +inf); div
  truncates. There is no per-call rounding-mode parameter.
- **Hann only.** No Hamming/Blackman tables yet; the `window()` function
  accepts any coefficient table, so adding one is a table, not new code.

## Declaration

- **Pre-existing work / reused code:** the ingot skeleton (`fe.toml`
  shape, `.gitignore`, CI structure) and the Apache-2.0 `LICENSE` are
  adapted from the official
  [builders-challenge-starter](https://github.com/fe-lang/builders-challenge-starter)
  (Apache-2.0). All DSP code in `src/lib.fe` and all tests in
  `src/tests.fe` are original, written for this challenge.
- **AI assistance:** this library was written with AI assistance (Muse).
  The author reviewed every line, verified the toolchain and language
  semantics by compiler probe, and verified every test vector against the
  independent Python oracle in `tools/`.
- **No prior publication:** this code has not been published or submitted
  anywhere before this challenge.

## License

Apache-2.0. See `LICENSE`.
