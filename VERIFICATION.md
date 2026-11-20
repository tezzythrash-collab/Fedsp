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
