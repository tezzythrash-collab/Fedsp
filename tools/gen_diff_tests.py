#!/usr/bin/env python3
"""Generate the v0.4.0 differential test functions from tools/vectors.json.

Mechanical transcription only: every expected value comes from the oracle
(tools/gen_vectors.py, seed 20261120). Re-run after regenerating
vectors.json, append the output to ingots/fedsp/src/tests.fe, then run
`fe fmt` on the test file.
"""
import json

v = json.load(open("tools/vectors.json"))
d = v["differential"]

QTYPE = {15: "Q15", 7: "Q7", 31: "Q31", 0: "Q<0>"}


def q15_list(vals, indent="        "):
    return (",\n" + indent).join(f"Q15::from_raw(raw: {x})" for x in vals)


def i32_list(vals):
    return ", ".join(str(x) for x in vals)


out = []
out.append("// --- v0.4.0 differential tests ----------------------------------------")
out.append("// Seeded random inputs (seed 20261120, see tools/gen_vectors.py).")
out.append("// Expected values transcribed from tools/vectors.json, never from")
out.append("// Fe output.")

# --- dot ---------------------------------------------------------------------
out.append("")
out.append("#[test]")
out.append("fn dot_random() {")
out.append("    // Random 8- and 16-sample frames through the saturating dot.")
for i, case in enumerate(d["dot"]):
    n = len(case["a"])
    out.append(f"    let a{i}: [Q15; {n}] = [")
    out.append(f"        {q15_list(case['a'])}")
    out.append("    ]")
    out.append(f"    let b{i}: [Q15; {n}] = [")
    out.append(f"        {q15_list(case['b'])}")
    out.append("    ]")
    out.append(f"    assert!(dot(a: a{i}, b: b{i}).raw() == {case['out']})")
out.append("}")

# --- FIR ---------------------------------------------------------------------
fir = d["fir"]
out.append("")
out.append("#[test]")
out.append("fn fir_random() {")
out.append("    // Random 8 taps over a random 16-sample input.")
out.append("    var fir = Fir::new(")
out.append("        taps: [")
out.append(f"            {q15_list(fir['taps'], indent='            ')}")
out.append("        ],")
out.append("    )")
out.append("    let input: [Q15; 16] = [")
out.append(f"        {q15_list(fir['input'])}")
out.append("    ]")
out.append("    let out = run_fir(fir: mut fir, xs: input)")
out.append("    assert_q15_frame(")
out.append("        got: out,")
out.append(f"        want: [{i32_list(fir['output'])}],")
out.append("    )")
out.append("}")

# --- biquad ------------------------------------------------------------------
bq = d["biquad"]
out.append("")
out.append("#[test]")
out.append("fn biquad_random() {")
out.append("    // Random 16-sample input through the documented lowpass.")
out.append("    let input: [Q15; 16] = [")
out.append(f"        {q15_list(bq['input'])}")
out.append("    ]")
out.append("    var lp = Biquad::lowpass_01()")
out.append("    let out = run_biquad(bq: mut lp, xs: input)")
out.append("    assert_q15_frame(")
out.append("        got: out,")
out.append(f"        want: [{i32_list(bq['output'])}],")
out.append("    )")
out.append("}")

# --- FFT ---------------------------------------------------------------------
out.append("")
out.append("#[test]")
out.append("fn fft_random() {")
out.append("    // Seeded random complex frames; the oracle cross-checked")
out.append("    // each against an independent float DFT (<= 16 ulp).")
for n in (16, 32):
    f = d["fft"][f"n{n}"]
    out.append(f"    let x{n} = to_cplx2(")
    out.append("        re: [")
    out.append(f"            {q15_list(f['in_re'], indent='            ')}")
    out.append("        ],")
    out.append("        im: [")
    out.append(f"            {q15_list(f['in_im'], indent='            ')}")
    out.append("        ],")
    out.append("    )")
    out.append(f"    let spec{n} = fft{n}(x: x{n})")
    out.append("    assert_cplx_frame(")
    out.append(f"        got: spec{n},")
    out.append(f"        want_re: [{i32_list(f['out_re'])}],")
    out.append(f"        want_im: [{i32_list(f['out_im'])}],")
    out.append("    )")
    out.append(f"    let mag{n} = fft_magnitude(x: spec{n})")
    out.append("    assert_q15_frame(")
    out.append(f"        got: mag{n},")
    out.append(f"        want: [{i32_list(f['magnitude'])}],")
    out.append("    )")
    out.append(f"    assert_parseval_cplx(x: x{n}, spec: spec{n})")
    out.append(f"    let rt{n} = ifft{n}(x: spec{n})")
    out.append("    assert_cplx_frame(")
    out.append(f"        got: rt{n},")
    out.append(f"        want_re: [{i32_list(f['roundtrip_re'])}],")
    out.append(f"        want_im: [{i32_list(f['roundtrip_im'])}],")
    out.append("    )")
out.append("}")

# --- scalars -----------------------------------------------------------------
out.append("")
out.append("#[test]")
out.append("fn scalar_random() {")
out.append("    // Seeded random raws (seed 20261120): sqrt, reciprocal,")
out.append("    // and convert_q over 6 format pairs.")
for raw, s, r in zip(d["scalars"], d["sqrt"], d["reciprocal"]):
    out.append(f"    assert!(Q15::from_raw(raw: {raw}).sqrt().raw() == {s})")
    out.append(f"    assert!(Q15::from_raw(raw: {raw}).reciprocal().raw() == {r})")
for i, c in enumerate(d["convert"]):
    ftype, ttype = QTYPE[c["from"]], QTYPE[c["to"]]
    out.append(f"    let c{i}: {ttype} = convert_q(x: {ftype}::from_raw(raw: {c['raw']}))")
    out.append(f"    assert!(c{i}.raw() == {c['out']})")
out.append("}")

print("\n".join(out))
