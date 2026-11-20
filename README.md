# FeDSP — Fixed-Point DSP Math Library for Fe

**One-sentence pitch:** Saturating binary Q-format fixed-point arithmetic for
audio-rate DSP in pure Fe — Q7/Q15/Q31, biquad and FIR filters, format
conversion, Hann/Hamming/Blackman windows, frame statistics, Newton
sqrt/reciprocal, and a radix-2 FFT / inverse FFT at 16/32/64 points, with
every result verified against an independent Python oracle.

**Fe-specific highlight:** the library is built on Fe 26.5 features that do
not exist in Solidity: a generic `Q<const FRAC: i32>` type (one
implementation serves every Q-format), a manual `Copy` impl that makes
fixed-size `Q` arrays usable with repetition and indexing, `i64` widening
with `downcast_saturate()` for overflow-free saturating arithmetic, a generic
FFT core (`fft_core<const N: usize>`) driven by oracle-generated twiddle and
bit-reversal tables, and pure functions that run unchanged on both the
sonatina and native test backends. The `DspDemo` contract shows the library
used as a real on-chain dependency: it runs the full signal chain inside a
`recv` handler and returns the results as a `#[abi]` struct.

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
  neg, dot products, filter steps, frame statistics, sqrt, and reciprocal
  all saturate to `i32::MAX` / `i32::MIN` instead of wrapping or
  reverting. Division by zero saturates by sign (or `0` for `0 / 0`).
  `dot` accumulates in `i64` with a saturating add at every step (v0.2.0
  fix — it used to revert on 3+ maximum-magnitude elements).
- **Exact rounding, documented.** Multiplication widens to `i64`, adds half
  an ulp, shifts, and clamps — round-to-nearest, ties toward +inf. Division
  truncates toward zero. `to_int` is an arithmetic shift (floors negatives:
  `-1.5 -> -2`). Every rule is stated in the doc comments and pinned by a
  test.
- **Verified, not asserted.** Every expected value in `ingots/fedsp/src/tests.fe` was
  computed by `tools/gen_vectors.py`, an independent Python oracle using
  exact arbitrary-precision integers — not by copying Fe output. If Fe and
  the oracle ever disagree, the test fails and the library is wrong.

## Architecture

```
fedsp/                          (Fe workspace, 0.4.0)
  fe.toml                  workspace: members = ["ingots/*"]
  ingots/fedsp/
    fe.toml                ingot: fedsp 0.4.0 (Fe 26.5.0+)
    src/lib.fe             the whole library (single file, ~1540 lines)
      Q<const FRAC: i32>   conversions, saturating arithmetic, scaling,
                           shifts, clamp/min/max, abs/neg, sqrt/reciprocal
      Q7 / Q15 / Q31       type aliases
      HANN16 / 32 / 64     Hann window tables, Q15 raw (Python-generated)
      HAMMING16 / 32 / 64  Hamming window tables, Q15 raw
      BLACKMAN16 / 32 / 64 Blackman window tables, Q15 raw
      dot<const FRAC, N>   dot product, saturating i64 accumulation,
                           rounded + clamped (never reverts)
      window<const N>      elementwise Hann (or any) windowing of Q15 frames
      convert_q<FROM, TO>  Q-format conversion, rounded + saturating
      mean / rms / energy  Q15 frame statistics, saturating i64 sums
      isqrt_round          integer sqrt, rounded to nearest (Newton)
      Biquad               transposed direct-form-II IIR, saturating step
      Fir<const N>         tap delay-line FIR over the saturating dot
      Cplx                 complex Q15 (re + i*im), saturating arithmetic
      TWIDDLE16/32/64_RE/IM  FFT twiddle factors W_N^k, Q15 raw (oracle-made)
      BITREV16 / 32 / 64   bit-reversal permutations for the FFT
      fft16 / 32 / 64      radix-2 DIT FFT (unscaled forward transform)
      ifft16 / 32 / 64     inverse FFT (conjugate trick, scaled by 1/N)
      fft_magnitude        |X[k]| per bin, Q15
    src/tests.fe           31 tests, 378 assertions, all oracle-verified
                           (26 fixed-vector tests + 5 seeded differential
                           tests, seed 20261120)
  ingots/demo/
    fe.toml                ingot: demo 0.4.0, depends on fedsp by path
    src/lib.fe             DspDemo contract + contract test (1 test,
                           sonatina only — see Limitations)
      ChainResult          #[abi] struct: rms, energy, fft_peak,
                           fft_peak_bin
      run_chain            the full demo chain as a pure function
      DspDemo              contract: RunDemo (no args) + RunChain
                           (caller-supplied 16-sample frame)
  tools/gen_vectors.py     the Python oracle (vectors + all window tables
                           + FFT twiddles, bit-reversal, demo chain,
                           differential vectors, demo-contract frame)
  tools/gen_diff_tests.py  transcribes the differential vectors into the
                           five differential test functions in
                           src/tests.fe (mechanical transcription only)
  tools/vectors.json       oracle output (all expected values)
  tools/window_tables.txt  oracle output (Hamming/Blackman tables, readable)
  tools/fft_tables.txt     oracle output (twiddle + bit-reversal tables)
  .github/workflows/ci.yml  fmt check, fe check, fe test (sonatina + native),
                           fe build
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
- Stateful filters use Fe 26.5's mutable-borrow pattern: `step(self: mut
  Self, ...)` mutates the state in place, the caller holds the filter in a
  `var` binding and passes `mut filter` to helpers (`fn run(bq: mut
  Biquad, ...)`). Plain parameters are read-only views — `var x = p`
  cannot move a parameter out, and `let mut` is not Fe (`var` is the
  mutable binding). `if`/`else` works as an expression; `usize as i64` is
  not allowed — use `.downcast_truncate()`.
- Message fields use braces (`RunChain { frame: [i32; 16] }`), and a
  struct used as a message return type needs `#[abi]` with `pub` fields
  (plain structs do not implement `Decode<Sol>`). Both were established
  by compiler probe while building the demo contract.
- A contract test poisons the native backend for its whole ingot:
  Fe 26.5.0's native backend cannot lower EVM instructions
  (`EvmSstore`, `EvmCallValue`), so `fe test --backend native` fails at
  codegen for any ingot containing a contract test — including unrelated
  pure tests in the same ingot (verified by probe). That is why the
  demo contract lives in its own workspace member (see Limitations).

## Where to look

Pointers to the evidence for each judging criterion — no sales copy,
just locations.

| Criterion | Weight | Where the evidence lives |
| --- | --- | --- |
| Technical quality and correctness | 30% | `ingots/fedsp/src/tests.fe` — 31 tests, 378 assertions, both backends; `tools/gen_vectors.py` — the independent oracle (fixed vectors + seeded differential vectors, seed 20261120); the `dot` overflow bug, found by probe and fixed with a regression test (`dot_saturation`), in `VERIFICATION.md` §v0.2.0; benchmarks in `VERIFICATION.md` §Benchmarks |
| Fe-specific value | 25% | `ingots/fedsp/src/lib.fe` — const-generic `Q<const FRAC: i32>` with compile-time bounds, manual `Copy` impl, `downcast_saturate()`, generic `fft_core<const N: usize>`, mutable-borrow stateful filters; `ingots/demo/src/lib.fe` — the `DspDemo` contract (`msg` + `recv` + `#[abi]` struct return) using the library as a dependency; the "Key design facts" list above (all compiler-probed Fe 26.5 behaviors) |
| Originality | 20% | `## Problem` above — the stdlib is decimal fixed-point for money; no binary Q-format DSP library existed for Fe. The combination (saturating Q arithmetic + windows + filters + FFT + oracle verification) is the original contribution |
| Usefulness or insight | 15% | `## Beyond audio: fixed-point math for DeFi` — sqrt/reciprocal as AMM primitives, Q-format vs decimal Wad/Ray; `ingots/demo/` — the library reused as a contract dependency by path (`fedsp = "../fedsp"`), the same way any future project would consume it |
| Documentation & demo | 10% | `## Demo` — the end-to-end signal chain with oracle-pinned values, plus the `DspDemo` contract (`run_demo()` / `run_chain(int32[16])`, both `pure` in `out/DspDemo.abi.json`); `VERIFICATION.md` — per-version verification records and benchmarks; `## Limitations` — stated plainly, including what is unmeasured |

## Beyond audio: fixed-point math for DeFi

FeDSP is written for audio-rate DSP, but its primitives are the same ones
on-chain math builds on — and that is deliberate, not an accident:

- **`sqrt` and `reciprocal` are AMM primitives.** Concentrated-liquidity
  math is built on sqrt-prices, and reciprocal is the other half of any
  price inversion. FeDSP's Newton `sqrt` (rounded to nearest, exact on the
  tested cases) and 5-step Newton `reciprocal` (~1–2 ulp) give contracts
  those operations in binary Q-format, next to the stdlib's decimal
  `Wad` / `Ray`.
- **Q-format complements decimal fixed-point.** The stdlib's `Fixed<D>` is
  decimal (`raw / 10^D`) — right for money. Binary Q-format (`raw /
  2^FRAC`) is right wherever a power of two maps to a shift: signal
  processing, sensor fusion, and any contract that ingests Q15 samples
  from hardware or codecs.
- **Saturation is the on-chain discipline.** A view or pure function that
  reverts on overflow is a denial-of-service vector for its callers;
  FeDSP saturates instead of reverting, everywhere, including the FFT
  butterflies — a filter or transform can never panic on adversarial
  input.

Honest limits, stated plainly: **gas use is unmeasured** — `fe test`
exercises the library off-chain, and no gas benchmark exists yet, so
nothing here is a claim about on-chain cost. FeDSP is a **math
primitive, not a DeFi product**: it ships no AMM, no vault, no oracle,
and no audit. The DeFi relevance is that the primitives are correct,
oracle-verified, and reusable — the product built on them is someone
else's (or a later version's) work.

## What's new in v0.4.0

- **Randomized differential tests.** Five new test functions
  (`dot_random`, `fir_random`, `biquad_random`, `fft_random`,
  `scalar_random`) run seeded random inputs (seed `20261120`, consumed
  in a fixed order so the vectors regenerate byte-identically) through
  `dot`, FIR, biquad, FFT (16/32), `sqrt`, `reciprocal`, and `convert_q`
  over 6 format pairs — 98 new assertions, on top of the 26 fixed-vector
  tests, not instead of them. The FFT random frames are cross-checked
  against an independent float DFT in the oracle (max error 5 ulp at 32
  points) with round-trips within 1 ulp. Expected values are transcribed
  mechanically from the oracle by `tools/gen_diff_tests.py`.
- **Benchmarks.** `VERIFICATION.md` now has a Benchmarks section:
  per-test wall-clock times on both backends from actual `fe test` runs
  on this VM, a complexity table for the main operations, and the
  `DspDemo` build artifact sizes. Framed honestly: wall-clock timings,
  not gas — gas use remains unmeasured.
- **Judges' map.** `## Where to look` above maps each judging criterion
  to the exact file and section where its evidence lives.
- **DspDemo contract.** The repo is now a workspace: the library lives
  in `ingots/fedsp/`, and `ingots/demo/` holds a small contract that
  runs the full signal chain on-chain — `run_demo()` (no arguments) and
  `run_chain(int32[16])` (caller-supplied frame), both `pure`, returning
  a `#[abi]` struct with RMS, energy, FFT peak magnitude, and peak bin.
  Its test deploys the contract in the EVM sandbox and asserts the
  oracle-pinned results for the demo burst and a seeded random frame.
  `fe build` now succeeds and produces `out/DspDemo.{bin,runtime.bin}`
  and `out/DspDemo.abi.json`. The demo contract's test is sonatina-only:
  Fe 26.5.0's native backend cannot lower EVM instructions (see
  Limitations) — which is why it lives in its own workspace member, so
  the library suite stays green on both backends.

## What's new in v0.3.0

- **Radix-2 FFT and inverse FFT.** `Cplx` (complex Q15, saturating
  arithmetic), `fft16` / `fft32` / `fft64` (iterative decimation-in-time,
  unscaled forward transform), `ifft16` / `ifft32` / `ifft64` (conjugate
  trick, scaled by exactly `1 / N`), and `fft_magnitude` (`|X[k]|` per
  bin). Twiddle factors and bit-reversal permutations are oracle-generated
  tables (`TWIDDLE16/32/64_RE/IM`, `BITREV16/32/64`). The oracle
  cross-checks the butterfly against an independent floating-point DFT
  (max error 11 ulp per component at 64 points) and pins
  forward-then-inverse round-trips to within 1 ulp, plus a Parseval check
  (`sum |X|^2 == N * sum |x|^2`, measured relative error ~2e-5).
- **End-to-end demo.** `demo_signal_chain` runs one signal through the
  whole library — Hann window, moving-average FIR, Butterworth lowpass
  biquad, RMS/energy, FFT magnitude — with every intermediate value
  pinned by the oracle (see Demo below). The old `demo_gain_stage` test
  still passes unchanged.
- **DeFi note.** The section above: what this library offers on-chain
  math, and what it does not claim.

## What's new in v0.2.0

- **`dot` overflow fix (the bug).** `dot` accumulated products with a plain
  `i64` add, which reverted with arithmetic-overflow panic (0x11) for 3+
  maximum-magnitude elements — breaking the library's headline rule
  ("saturate, never revert"). The accumulator now uses
  `acc.saturating_add(product)` at every step, and the final half-add for
  rounding is saturating too. Regression tests pin `dot` of 3 and 8
  maximum-magnitude elements (both signs) to `maximum()` / `minimum()`.
- **Biquad filter.** Transposed direct-form II, Q15 coefficients
  (b0/b1/b2/a1/a2) and state (z1/z2), saturating step, `reset()`, and a
  documented `lowpass_01()` example (2nd-order Butterworth, fc = 0.1 fs;
  coefficients [2210, 4421, 2210, -37453, 13527]). Impulse and
  step responses are oracle-verified.
- **FIR filter.** Generic tap delay line over the fixed `dot`, with
  oracle-verified impulse/step responses (moving average + random 8 taps).
- **Format conversion.** `convert_q<FROM, TO>` between any Q-formats
  (including `Q<0>`), rounding to nearest (ties toward +inf) and
  saturating; tested across Q7/Q15/Q31 boundaries.
- **Hamming and Blackman tables.** `HAMMING16/32/64` and
  `BLACKMAN16/32/64`, generated by the oracle in the same style as the
  Hann tables, with sum/peak/endpoint tests.
- **Frame statistics.** `mean`, `rms`, and `energy` for Q15 frames,
  saturating `i64` accumulation, oracle-verified (including the
  ties-toward-+inf rounding: mean of `[-1, -2]` is `-1`, not `-2`).
- **Newton sqrt and reciprocal.** `Q::sqrt` (integer Newton sqrt, rounded
  to nearest; negative input returns zero) and `Q::reciprocal` (5 fixed
  Newton steps on a normalized input, ~1-2 ulp of the true value; `0`
  returns `maximum()`; `Q<31>` saturates by sign since every reciprocal
  there is unrepresentable; `Q<0>` gives the truncated integer
  reciprocal). Edge cases are tested on both backends.

## Demo (reproducible CLI example)

The `demo_signal_chain` test is the demo: a 16-sample sine burst (2 cycles,
amplitude 0.9) goes through a Hann window, a 4-tap moving-average FIR, a
Butterworth lowpass biquad (fc = 0.1 fs), an RMS/energy measurement, and an
FFT magnitude — the whole library in one pipeline, in pure Fe, with
oracle-verified values at each step (filtered RMS `5477`, energy `14649`,
FFT magnitude peaks `42578` / `42686` at bins 1 and 2, mirror image at bins
14 and 15). The earlier `demo_gain_stage` test (window, 0.5 gain, energy)
still passes unchanged.

### The demo as a contract

The same chain also runs on-chain. `ingots/demo/` holds the `DspDemo`
contract, which depends on the library by path and exposes two `pure`
functions (see `out/DspDemo.abi.json` after `fe build`):

- `run_demo()` — runs the chain on the documented demo burst, no
  arguments. Returns `(rms, energy, fft_peak, fft_peak_bin)` =
  `(5477, 14649, 42686, 2)`.
- `run_chain(int32[16])` — runs the chain on a caller-supplied 16-sample
  Q15 frame. The contract test calls it with the demo burst (same result
  as `run_demo`) and with a seeded random frame (seed 20261120), where it
  must return `(552, 149, 4934, 2)` — the oracle-pinned values for that
  frame.

The contract test deploys `DspDemo` in the EVM sandbox and asserts all
three calls. It runs on the sonatina backend only (see Limitations).

```console
$ fe test
PASS  [0.0000s]  fedsp conversions_from_i64
PASS  [0.0000s]  fedsp conversions_from_int
PASS  [0.0000s]  fedsp conversions_from_ratio
PASS  [0.0000s]  fedsp convert_q_formats
PASS  [0.0000s]  fedsp q0_plain_integers
PASS  [0.0000s]  fedsp sat_add_sub
PASS  [0.0000s]  fedsp sat_mul
PASS  [0.0001s]  fedsp abs_neg_clamp_min_max
PASS  [0.0001s]  fedsp constants
PASS  [0.0001s]  fedsp cplx_arithmetic
PASS  [0.0001s]  fedsp demo_gain_stage
PASS  [0.0001s]  fedsp dot_products
PASS  [0.0001s]  fedsp dot_random
PASS  [0.0001s]  fedsp dot_saturation
PASS  [0.0001s]  fedsp sat_div
PASS  [0.0001s]  fedsp scale_and_shifts
PASS  [0.0002s]  fedsp biquad_random
PASS  [0.0002s]  fedsp fir_random
PASS  [0.0003s]  fedsp frame_stats
PASS  [0.0003s]  fedsp hann_tables
PASS  [0.0004s]  fedsp biquad_filter
PASS  [0.0005s]  fedsp scalar_random
PASS  [0.0005s]  fedsp window_tables_hamming_blackman
PASS  [0.0007s]  fedsp sqrt_reciprocal
PASS  [0.0013s]  fedsp demo_signal_chain
PASS  [0.0017s]  fedsp fft16_transform
PASS  [0.0020s]  fedsp fir_filter
PASS  [0.0048s]  fedsp fft_random
PASS  [0.0053s]  fedsp fft32_transform
PASS  [0.0059s]  fedsp fft64_transform
PASS  [0.0060s]  fedsp ifft_roundtrip
PASS  [0.0097s]  demo  test_dsp_demo

test result: ok. 32 passed; 0 failed
```

(Timings are from one run; they vary. The pass list and the result line
are what matter.) The 31 library tests also pass on the native backend
(`fe test --backend native ingots/fedsp`); the demo contract test is
sonatina-only (see Limitations).

## Build & test

Requires Fe **26.5.0** or later (verified with `fe 26.5.0 (3dab95f)`).

```console
# Install the toolchain (Linux x86_64; also works on macOS/Windows via feup)
curl -fsSL https://raw.githubusercontent.com/argotorg/fe/master/feup/feup.sh -o /tmp/feup.sh
bash /tmp/feup.sh --version v26.5.0 --no-modify-path
export PATH="$HOME/.fe/bin:$PATH"

# In the repo root (a workspace; members are ingots/fedsp + ingots/demo):
fe fmt --check ingots/fedsp   # formatting, library
fe fmt --check ingots/demo    # formatting, demo contract
fe check                      # type check, whole workspace
fe test                       # 32 tests, sonatina backend (both members)
fe test --backend native ingots/fedsp   # 31 library tests, native backend
fe build                      # compiles the DspDemo contract to out/
```

`fe build` produces `out/DspDemo.bin` (deploy bytecode), `out/DspDemo.runtime.bin`,
and `out/DspDemo.abi.json`. `out/` is git-ignored; the artifacts are build
outputs, not committed sources.

Note: `fe test --backend native` runs the library ingot only. The demo
contract's test is sonatina-only — Fe 26.5.0's native backend cannot lower
EVM instructions (see Limitations). `fe test` at the workspace root runs
both members on sonatina.

To regenerate the test vectors and window tables:

```console
python3 tools/gen_vectors.py   # rewrites tools/vectors.json, tools/hann_tables.txt, tools/window_tables.txt, tools/fft_tables.txt
```

## Limitations

- **i32 raw only.** No `i64`-raw or `u32` Q-formats yet; `Q31` tops out
  just below `1.0` (see `one()` above).
- **FRAC range `0..=31`, enforced at compile time.** `Q` carries a
  `where FRAC >= 0, FRAC <= 31` bound (propagated to `dot`); instantiating
  `Q` outside the range is a compile error, not undefined behavior.
- **FFT sizes are 16, 32, and 64 only.** The transform is radix-2, so
  only power-of-two sizes work, and only the three sizes with generated
  twiddle tables are exposed. Other sizes need new tables (the oracle
  generates them) and new wrappers.
- **The forward FFT is unscaled.** A bin magnitude can grow by up to `N`
  over the input magnitude. For real audio inputs (`|x| <= 1.0` in Q15)
  the largest bin is `N * 32768` raw, far inside `i32`; saturation only
  bites for near-full-range `i32` raws, and then the butterflies clamp
  instead of reverting. The inverse FFT scales by exactly `1 / N`, so
  `ifft(fft(x)) == x` within 1 ulp. The Q15 twiddles themselves are
  rounded to the nearest raw, which is why the transform agrees with an
  exact float DFT to ~11 ulp per component at 64 points, not exactly.
- **`fft_magnitude` saturates above |X| = 256.** The magnitude is exact
  while a bin magnitude is at most 256 (Q15 value); above that the
  squared sum saturates. For real audio inputs the largest bin is `N`
  (64 at most), so this never triggers at 16/32/64 points.
- **No resampling (yet).** The library now covers the arithmetic core,
  windowing, dot products, biquad/FIR filters, format conversion, frame
  statistics, sqrt, reciprocal, and the FFT — resampling is the remaining
  natural follow-up, not included here.
- **Reciprocal is ~1-2 ulp, not exact.** Newton's iteration converges to
  within ~1-2 ulp of the true reciprocal (the rounded iteration stalls 1
  ulp below values like 2.0 and 4.0); the exact result is pinned by the
  oracle, but callers needing a correctly-rounded reciprocal should not
  assume it. For small inputs the absolute error grows with the final
  left shift, while the relative error stays ~1 ulp.
- **`Q<31>` reciprocal always saturates.** Every `|1/x| > 1` for
  `|x| < 1` in Q31, so `reciprocal` returns `maximum()` / `minimum()` by
  sign (except `reciprocal(0) == maximum()`).
- **Frame statistics are Q15-only.** `mean`, `rms`, and `energy` take
  `[Q15; N]`; other formats need `convert_q` first.
- **The demo contract's test is sonatina-only.** Fe 26.5.0's native
  backend cannot lower EVM instructions (`EvmSstore`, `EvmCallValue`),
  so any ingot containing a contract test fails `fe test --backend
  native` at codegen — including unrelated pure tests in the same
  ingot (verified by probe: a stateless contract test and a pure test
  in separate files still failed native together, 0 passed). That is
  why `DspDemo` lives in its own workspace member (`ingots/demo/`):
  the library suite stays green on both backends, and the contract
  test runs on sonatina via `fe test ingots/demo` (or `fe test` at the
  workspace root). If a future Fe release lowers EVM instructions on
  the native backend, the demo test should be run there too.
- **No `f32` interop.** Conversion to/from float is out of scope for this
  version; all boundaries are integer.
- **Rounding is fixed.** Mul rounds to nearest (ties toward +inf); div
  truncates. There is no per-call rounding-mode parameter.

## Declaration

- **Pre-existing work / reused code:** the ingot skeleton (`fe.toml`
  shape, `.gitignore`, CI structure) and the Apache-2.0 `LICENSE` are
  adapted from the official
  [builders-challenge-starter](https://github.com/fe-lang/builders-challenge-starter)
  (Apache-2.0). All DSP code in `ingots/fedsp/src/lib.fe`, all tests in
  `ingots/fedsp/src/tests.fe`, and the `DspDemo` contract in
  `ingots/demo/src/lib.fe` are original, written for this challenge.
  The workspace layout and the demo contract were added in v0.4.0,
  inside the contest window.
- **AI assistance:** this library was written with AI assistance (Muse).
  The author reviewed every line, verified the toolchain and language
  semantics by compiler probe, and verified every test vector against the
  independent Python oracle in `tools/`.
- **No prior publication:** this code has not been published or submitted
  anywhere before this challenge.

## License

Apache-2.0. See `LICENSE`.
