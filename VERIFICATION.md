# Verification — FeDSP v0.4.0 (2026-11-20)

## What changed

- **Workspace restructure.** The repo is now a Fe workspace: the library
  ingot moved to `ingots/fedsp/` (renames only, `git mv`), and
  `ingots/demo/` holds the DspDemo contract ingot. Root `fe.toml` is the
  workspace config (`members = ["ingots/*"]`). `fe.toml` versions bumped
  to 0.4.0 (workspace + both members).
- **Randomized differential tests** (seed `20261120`, consumed in the
  fixed order dot -> fir -> biquad -> fft -> scalars -> demo-contract
  frame, so regeneration is byte-identical). Five new test fns
  (`dot_random`, `fir_random`, `biquad_random`, `fft_random`,
  `scalar_random`), 98 new `assert!` sites: random 8/16-sample frames
  through `dot` (6 pairs), random 8 taps + random input through FIR,
  random input through the lowpass biquad, random complex frames through
  FFT 16/32 (oracle cross-check vs independent float DFT: max error 2
  ulp at 16 points, 5 ulp at 32 — bound 16; round-trips within 1 ulp —
  bound 2), 16 random raws through `sqrt`/`reciprocal`, and `convert_q`
  over 6 format pairs for 10 random raws (60 assertions). Expected values
  are transcribed mechanically from the oracle by
  `tools/gen_diff_tests.py` — the script reads `tools/vectors.json` and
  emits the test functions; nothing was copied from Fe output.
- **DspDemo contract** (`ingots/demo/src/lib.fe`): `run_demo()` and
  `run_chain(int32[16])`, both `pure`, returning a `#[abi]` struct
  (`rms`, `energy`, `fft_peak`, `fft_peak_bin`). The test deploys the
  contract in the EVM sandbox and asserts the oracle-pinned results for
  the demo burst `(5477, 14649, 42686, 2)` and a seeded random frame
  `(552, 149, 4934, 2)`. `fe build` now succeeds:
  `out/DspDemo.bin` 12,940 B, `out/DspDemo.runtime.bin` 12,900 B,
  `out/DspDemo.abi.json` 1,201 B (out/ is git-ignored).
- **Judges' map** in README.md (`## Where to look`): each rubric
  criterion mapped to the exact file/section of its evidence.
- **Benchmarks** — see the Benchmarks section below.

## What was run (this VM, `fe 26.5.0 (3dab95f)`, in `~/workspace/fedsp/`)

- `fe fmt --check ingots/fedsp` → clean. `fe fmt --check ingots/demo` → clean.
- `fe check` (workspace root) → clean.
- `fe test` (sonatina, both members) → **32 passed, 0 failed**, exit 0.
- `fe test --backend native ingots/fedsp` → **31 passed, 0 failed**, exit 0.
- `fe build` (workspace root) → succeeds, writes the three DspDemo artifacts.
- New test names: dot_random, fir_random, biquad_random, fft_random,
  scalar_random (library), test_dsp_demo (demo contract). Total: 32 test
  fns, 391 `assert!` sites (library 378, demo contract 13).
- Oracle: `python3 tools/gen_vectors.py` re-run; all v0.3.0 vectors
  byte-identical (verified by diff — additions only: `differential`
  and `demo_contract` keys); all four table files unchanged. New oracle
  sections: differential vectors (seed 20261120) and the demo-contract
  frame. New tool: `tools/gen_diff_tests.py`.

## Native-backend limitation (verified by probe, not assumed)

Fe 26.5.0's native backend cannot lower EVM instructions. Probes in
`/tmp/fe-probe` (scratch ingots, repo untouched):
- The `fe new` Counter scaffold: `fe test` passes on sonatina;
  `fe test --backend native` fails codegen on `EvmSstore` in
  `__Counter_init`.
- A stateless contract (no storage, pure `recv` handler): native fails
  codegen on `EvmCallValue` in `contract_init_abi_Demo`.
- A stateless contract test plus a pure test in a separate file: native
  fails the whole run at `tests::codegen` (0 passed) — the pure test
  never runs. Sonatina passes both.
- The workspace probe: `fe test ingots/lib` passes on both backends,
  `fe test ingots/demo` passes on sonatina, `fe build` writes the
  artifacts, `fe test` at the workspace root runs all members.
- Message-syntax probes: msg variant fields use braces
  (`RunChain { frame: [i32; 16] }`); a struct return type needs `#[abi]`
  with `pub` fields (a plain struct fails with "`ChainResult` doesn't
  implement `Decode<Sol>`").

Consequence, applied in this repo: the demo contract lives in its own
workspace member, the library suite stays green on both backends, and
the contract test runs on sonatina only. This is a toolchain limitation
of Fe 26.5.0, not a property of the library or the contract.

## Audit notes

- `ingots/fedsp/src/lib.fe` (1537 lines, unchanged by v0.4.0): pure
  functions + two stateful filter structs + the FFT, no storage, no EVM
  effects, no network, no new dependencies.
- `ingots/fedsp/src/tests.fe`: every expected value comes from the
  oracle; the differential tests add randomized coverage on top of the
  fixed vectors, transcribed mechanically by `tools/gen_diff_tests.py`.
- `ingots/demo/src/lib.fe`: pure computation only (no storage, no
  effects beyond the call); the contract imports the library by path
  (`fedsp = "../fedsp"`), the same way any future project would consume
  it. Both messages are `pure` in the ABI.
- Known limits (also in README.md): FFT sizes 16/32/64 only; forward FFT
  unscaled (documented saturation behavior); `fft_magnitude` saturates
  above |X| = 256 (never triggers for real audio inputs at these sizes);
  no resampling yet; gas use unmeasured (see Benchmarks).
- Status: **VERIFIED** (library both backends; demo contract sonatina).
  **NOT pushed, NOT re-submitted.** The FE-089 form still points at the
  v0.1.0 head commit `0a0159006392f1e14cd4351fa12653c57d63c817`;
  re-submission before Nov 24, 2026, 23:59 UTC is a separate founder step
  if v0.4.0 is to be judged.
- Post-verification fix (independent review): README.md "Verified, not
  asserted" cited the pre-workspace path `src/tests.fe`; corrected to
  `ingots/fedsp/src/tests.fe`. Docs-only change; `fe test` re-run after
  the fix: 32/32 sonatina.
- Code commit: `8353a956` (README path fix; this VERIFICATION.md committed separately as
  the head; the two commits differ only in this file). The v0.4.0 code commit was
  `40cbbe14850da511cfc7140b0b3262afa1cb27da`.

## Benchmarks

Wall-clock timings from actual `fe test` runs on this VM (single runs;
times vary run to run, and the native backend runs each test in its own
process, so its numbers include per-test startup overhead). These are
**wall-clock timings, not gas**: gas use remains unmeasured — `fe test`
exercises the library off-chain, and no gas benchmark exists yet.

Per-test wall-clock times, sonatina backend (31 library tests, total
0.0289 s; slowest: `fft64_transform` 0.0082 s, `ifft_roundtrip`
0.0057 s, `fft_random` 0.0049 s; every other test ≤ 0.0026 s):

| Test | Time |
| --- | --- |
| fft64_transform | 0.0082 s |
| ifft_roundtrip | 0.0057 s |
| fft_random | 0.0049 s |
| fft32_transform | 0.0026 s |
| fft16_transform | 0.0013 s |
| demo_signal_chain | 0.0012 s |
| sqrt_reciprocal | 0.0007 s |
| biquad_filter | 0.0006 s |
| scalar_random | 0.0006 s |
| window_tables_hamming_blackman | 0.0005 s |
| fir_filter | 0.0004 s |
| frame_stats | 0.0003 s |
| biquad_random | 0.0002 s |
| fir_random | 0.0002 s |
| demo_gain_stage | 0.0002 s |
| hann_tables | 0.0002 s |
| all other library tests (15) | ≤ 0.0001 s each |

Per-test wall-clock times, native backend (31 library tests, total
1.0788 s; slowest: `convert_q_formats` 0.0892 s, `demo_gain_stage`
0.0749 s, `window_tables_hamming_blackman` 0.0744 s — the spread across
tests is dominated by per-test process overhead, not by the operations
under test, so the native column is a smoke signal, not a ranking):

| Test | Time |
| --- | --- |
| convert_q_formats | 0.0892 s |
| demo_gain_stage | 0.0749 s |
| window_tables_hamming_blackman | 0.0744 s |
| biquad_filter | 0.0740 s |
| frame_stats | 0.0692 s |
| dot_saturation | 0.0640 s |
| q0_plain_integers | 0.0630 s |
| cplx_arithmetic | 0.0563 s |
| fft16_transform | 0.0461 s |
| constants | 0.0418 s |
| sat_div | 0.0363 s |
| hann_tables | 0.0350 s |
| abs_neg_clamp_min_max | 0.0295 s |
| scale_and_shifts | 0.0291 s |
| sat_mul | 0.0289 s |
| fir_random | 0.0269 s |
| sat_add_sub | 0.0258 s |
| fft32_transform | 0.0255 s |
| sqrt_reciprocal | 0.0252 s |
| fft_random | 0.0252 s |
| ifft_roundtrip | 0.0216 s |
| scalar_random | 0.0208 s |
| dot_random | 0.0201 s |
| demo_signal_chain | 0.0165 s |
| conversions_from_ratio | 0.0109 s |
| conversions_from_i64 | 0.0098 s |
| fft64_transform | 0.0098 s |
| dot_products | 0.0087 s |
| fir_filter | 0.0079 s |
| biquad_random | 0.0065 s |
| conversions_from_int | 0.0059 s |

Demo contract test (sonatina): `test_dsp_demo` 0.0037 s (0.0097 s on
the workspace-root run — varies). `fe build` writes `DspDemo.bin`
12,940 B, `DspDemo.runtime.bin` 12,900 B, `DspDemo.abi.json` 1,201 B.

Complexity (per operation, from the implementation in
`ingots/fedsp/src/lib.fe` — operation counts, not measured times):

| Operation | Complexity | Notes |
| --- | --- | --- |
| `sat_add` / `sat_sub` / `sat_mul` / `sat_div` | O(1) | one i64 widen + round + clamp |
| `convert_q` | O(1) | one shift + rounding add |
| `Q::sqrt` | O(log bits) Newton steps | `isqrt_round`: Newton from `x/2 + 1`, converges while the next iterate is smaller; bounded by ~31 iterations for i64 |
| `Q::reciprocal` | O(1), 5 fixed Newton steps | plus a normalization loop of at most 31 shifts |
| `dot<const N>` | O(N) | N `sat_mul` + N saturating i64 adds + 1 rounding add |
| `window` | O(N) | N `sat_mul` |
| `mean` / `energy` | O(N) | N saturating i64 adds; `energy` is `dot(x, x)` |
| `rms` | O(N) + 1 `isqrt_round` | saturating sum of squares, then one integer sqrt |
| `Biquad::step` | O(1) per sample | 5 `sat_mul` + 3 saturating add/sub |
| `Fir::step` | O(N) per sample | one `dot` over N taps |
| `fft_core<const N>` | O(N log N) butterflies | each butterfly: 1 `Cplx::mul` (4 `sat_mul`) + 2 `Cplx` add/sub; N = 16/32/64 |
| `fft_magnitude` | O(N) | one magnitude (2 squares + sqrt) per bin |

---

# Verification — FeDSP (Batch A re-run, 2026-11-18)

## What was run (independent re-run, this VM)

- `fe 26.5.0 (3dab95f)` at `~/.fe/bin/fe`, in `~/workspace/fedsp/`.
- `fe test` (default sonatina backend) → **13 passed, 0 failed**, exit 0.
- `fe test --backend native` → **13 passed, 0 failed**, exit 0.
- Test names: conversions_from_int, conversions_from_i64, conversions_from_ratio,
  constants, abs_neg_clamp_min_max, q0_plain_integers, sat_add_sub, sat_mul,
  sat_div, scale_and_shifts, dot_products, demo_gain_stage, hann_tables.
- Submission package `~/workspace/your_files/fedsp-submission-package/fedsp-submission-package.zip`
  (21,987 B, 17 files) listed only — contents not extracted, nothing submitted.
  The private submission link/code was not touched, recorded, or pasted.

## Audit notes

- `src/lib.fe`: pure functions, no storage, no EVM effects; saturating arithmetic
  throughout; division by zero saturates by sign (0 for 0/0) instead of reverting.
- `src/tests.fe`: 13 test fns, expected values computed by the independent Python
  oracle `tools/gen_vectors.py` (exact arbitrary-precision integers).
- No network, no spend, no stubs, no TODO/FIXME markers found in `src/`, `tools/`.
- Status: **VERIFIED** (both backends). NOT submitted — submission remains a
  founder step; deadline on record is Nov 24, 2026, 23:59 UTC.

---

# Verification — FeDSP v0.2.0 (2026-11-20)

## What changed

- **`dot` overflow bug fixed.** `dot` accumulated products with a plain
  `i64` add and reverted with arithmetic-overflow panic (0x11) for 3+
  maximum-magnitude elements — proven by probe on both backends before the
  fix. The accumulator now uses `acc.saturating_add(product)` at every
  step, and the final half-add for rounding is saturating too. The `dot`
  doc comment no longer claims "no intermediate overflow for realistic
  lengths"; it states the saturation behavior. Regression tests
  (`dot_saturation`) pin `dot` of 3 and 8 maximum-magnitude elements, both
  signs and mixed signs, to `maximum()` / `minimum()`.
- **Biquad filter** (transposed direct form II, Q15, saturating step,
  `lowpass_01()` example), **FIR filter** (tap delay line over the fixed
  `dot`), **`convert_q`** (Q-format conversion, rounded + saturating),
  **Hamming/Blackman tables** (16/32/64), **`mean`/`rms`/`energy`** (Q15
  frame statistics), **`sqrt`/`reciprocal`** (Newton iteration).
- `Biquad::step` and `Fir::step` use Fe 26.5's mutable-borrow pattern
  (`step(self: mut Self, ...)`); the caller holds the filter in a `var`
  binding. This was established by compiler probe, not assumed: plain
  parameters are read-only views and cannot be moved into a `var`.

## What was run (this VM, `fe 26.5.0 (3dab95f)`, in `~/workspace/fedsp/`)

- `fe fmt --check` → clean. `fe check` → clean.
- `fe test` (sonatina) → **20 passed, 0 failed**, exit 0.
- `fe test --backend native` → **20 passed, 0 failed**, exit 0.
- New test names: dot_saturation, convert_q_formats,
  window_tables_hamming_blackman, frame_stats, sqrt_reciprocal,
  biquad_filter, fir_filter. Total: 20 test fns, 222 `assert!` sites
  (v0.1.0: 13 fns, 105 sites).
- Oracle: `python3 tools/gen_vectors.py` re-run; `tools/vectors.json`
  diff is additions only (460 lines added, 0 deleted — all v0.1.0 vectors
  byte-identical); `tools/hann_tables.txt` unchanged. New oracle sections:
  saturating dot, convert_q, q_sqrt, q_reciprocal, mean/rms/energy,
  biquad (Butterworth lowpass fc=0.1fs, coefficients verified by hand:
  [2210, 4421, 2210, -37453, 13527]), FIR, Hamming/Blackman tables.
  New file `tools/window_tables.txt` (oracle output, human-readable).
- Two test bugs found and fixed during this work (both were wrong expected
  values in the tests, not library bugs): the gain-biquad test expected
  8192 for a 0.5 gain on 1.0 (correct: 16384), and the oracle's `mean`
  rounding rule for negative sums with even `n` rounded a -0.5 tie the
  wrong way (mean of `[-1, -2]` is `-1`, not `-2` — ties toward +inf, the
  same convention as `sat_mul`). The oracle rule and the Fe `mean` were
  both corrected to `rem <= -(n/2 + 1)`.

## Audit notes

- `src/lib.fe` (935 lines, was 427): pure functions + two stateful filter
  structs, no storage, no EVM effects, no network, no new dependencies.
- `src/tests.fe` (731 lines): every expected value comes from the oracle;
  nothing was copied from Fe output. The biquad/FIR impulse and step
  responses, the window-table sums/peaks, and the sqrt/reciprocal edge
  cases are all oracle-pinned.
- Known limits (also in README.md): reciprocal is ~1-2 ulp, not exact;
  `Q<31>` reciprocal always saturates; frame statistics are Q15-only;
  no FFT/resampling yet.
- Status: **VERIFIED** (both backends). **NOT pushed, NOT re-submitted.**
  The FE-089 form still points at the v0.1.0 head commit
  `0a0159006392f1e14cd4351fa12653c57d63c817`; re-submission before Nov 24,
  2026, 23:59 UTC is a separate founder step if v0.2.0 is to be judged.
- Code commit: `c9e310582ff7c91680173cb1fd61fe4a23539707` (this VERIFICATION.md committed
  separately as the head; the two commits differ only in this file).

---

# Verification — FeDSP v0.3.0 (2026-11-20)

## What changed

- **Radix-2 FFT and inverse FFT** (`Cplx`, `fft16` / `fft32` / `fft64`,
  `ifft16` / `ifft32` / `ifft64`, `fft_magnitude`). Iterative
  decimation-in-time over Q15 complex samples, driven by
  oracle-generated twiddle tables (`TWIDDLE16/32/64_RE/IM`) and
  bit-reversal permutations (`BITREV16/32/64`). The forward transform is
  unscaled (bins grow by up to `N`); the inverse scales by exactly `1/N`
  with a rounding divide (`div_pow2_round`, ties toward +inf). Every
  butterfly multiply/add saturates, so the transform never reverts.
- **End-to-end demo** (`demo_signal_chain`): sine burst -> Hann window ->
  4-tap moving-average FIR -> Butterworth lowpass biquad -> RMS/energy ->
  FFT magnitude, every intermediate value oracle-pinned. The old
  `demo_gain_stage` test is unchanged and still passes.
- **DeFi usefulness note** in README.md: sqrt/reciprocal as AMM
  primitives, Q-format vs decimal Wad/Ray, saturation as the on-chain
  discipline — with the explicit limits that gas use is unmeasured and
  the library is a math primitive, not a DeFi product.
- `fe.toml` bumped to 0.3.0.

## What was run (this VM, `fe 26.5.0 (3dab95f)`, in `~/workspace/fedsp/`)

- `fe fmt --check` → clean. `fe check` → clean.
- `fe test` (sonatina) → **26 passed, 0 failed**, exit 0.
- `fe test --backend native` → **26 passed, 0 failed**, exit 0.
- New test names: cplx_arithmetic, fft16_transform, fft32_transform,
  fft64_transform, ifft_roundtrip, demo_signal_chain. Total: 26 test fns,
  280 `assert!` sites (v0.2.0: 20 fns, 222 sites).
- Oracle: `python3 tools/gen_vectors.py` re-run; all v0.2.0 vectors
  byte-identical (verified by diff against the v0.2.0 output — additions
  only). New oracle sections: FFT twiddles + bit-reversal, the mirrored
  butterfly, the float-DFT cross-check, the inverse FFT, the demo chain,
  `Cplx` unit values. New file `tools/fft_tables.txt` (oracle output,
  human-readable).
- Oracle cross-checks (fail loudly at generation time if violated):
  mirrored butterfly vs independent float DFT, max error **11 ulp** per
  component at 64 points (bound 16); forward-then-inverse round-trip,
  max error **1 ulp** (bound 2); Parseval `sum |X|^2 == N * sum |x|^2`,
  measured relative error ~2e-5 (Fe test bound 1e-3, with wide headroom).
- One Fe type-inference fix during this work (not a library bug):
  `N.downcast_truncate() * tsum` in the Parseval helpers had to become
  two statements (`let n: i64 = N.downcast_truncate(); let expect = n *
  tsum`) — the same shape `mean()` in `src/lib.fe` already uses.

## Audit notes

- `src/lib.fe` (1537 lines, was 935): pure functions + two stateful
  filter structs + the FFT, no storage, no EVM effects, no network, no
  new dependencies. The formatter (`fe fmt`) expanded the new twiddle /
  bit-reversal tables and the new test arrays to one value per line; the
  pre-existing tables were not reformatted.
- `src/tests.fe` (2524 lines): every expected value comes from the oracle;
  nothing was copied from Fe output. The FFT impulse / DC / sine / random
  frames, the round-trip values, and the demo-chain intermediates are all
  oracle-pinned; the Parseval checks are structural (they would fail by
  far more than the 1e-3 bound for a transform that is not the DFT).
- Known limits (also in README.md): FFT sizes 16/32/64 only; forward FFT
  unscaled (documented saturation behavior); `fft_magnitude` saturates
  above |X| = 256 (never triggers for real audio inputs at these sizes);
  no resampling yet; gas use unmeasured.
- Nothing left out: all three requested upgrades (FFT, demo, DeFi note)
  are built and green. No improvement proved infeasible in Fe 26.5.
- Status: **VERIFIED** (both backends). **NOT pushed, NOT re-submitted.**
  The FE-089 form still points at the v0.1.0 head commit
  `0a0159006392f1e14cd4351fa12653c57d63c817`; re-submission before Nov 24,
  2026, 23:59 UTC is a separate founder step if v0.3.0 is to be judged.
- Code commit: `00da348146829f0c4f8973a8dc9ea53c8200a819` (this VERIFICATION.md committed
  separately as the head; the two commits differ only in this file).
