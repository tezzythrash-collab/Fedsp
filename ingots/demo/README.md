# DspDemo — the FeDSP library as an on-chain dependency

A small contract that runs the full FeDSP signal chain inside a `recv`
handler: Hann window → 4-tap moving-average FIR → Butterworth lowpass
biquad → RMS / energy → FFT magnitude peak. It exists to prove the
library works as a dependency of real contract code, not only in unit
tests.

## Interface

Both functions are `pure` (see `../../out/DspDemo.abi.json` after `fe build`):

- `run_demo()` — runs the chain on the documented demo burst (16
  samples, 2 cycles, amplitude 0.9), no arguments. Returns
  `(rms, energy, fft_peak, fft_peak_bin)` = `(5477, 14649, 42686, 2)`.
- `run_chain(int32[16])` — runs the chain on a caller-supplied
  16-sample Q15 frame. Returns the same `#[abi]` struct.

The contract test (`test_dsp_demo`, in `src/lib.fe`) deploys the
contract in the EVM sandbox and asserts all three calls: `run_demo`,
`run_chain` with the demo burst, and `run_chain` with a seeded random
frame (seed 20261120), which must return `(552, 149, 4934, 2)`. Every
expected value is pinned by the oracle (`../../tools/vectors.json`,
`demo_chain` / `demo_contract`).

## Commands (from the repo root)

```console
fe test ingots/demo    # the contract test, sonatina backend
fe build               # compiles DspDemo to out/
```

## Why this is a separate ingot

Fe 26.5.0's native backend cannot lower EVM instructions
(`EvmSstore`, `EvmCallValue`), so any ingot containing a contract test
fails `fe test --backend native` at codegen — including unrelated pure
tests in the same ingot (verified by probe). Keeping the contract in
its own workspace member lets the library suite stay green on both
backends while the contract test runs on sonatina. If a future Fe
release lowers EVM instructions on the native backend, the demo test
should be run there too.
