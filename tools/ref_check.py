#!/usr/bin/env python3
"""Differential check of ribbit-core against the C++ reference (Rattlegram), and golden export.

  ref_check.py --core build/ribbit_tool --ref build/ref_tool --work build/refcheck
      TX: every vector rendered by both, compared byte for byte.
      RX: clean, noisy (AWGN), frequency-shifted, level-scaled and truncated channels decoded
          by both; the event lists (kind, mode, call, CFO as hex float, polar flips, payload,
          sample index) must be identical. Clean channels must decode.
  ref_check.py --ref build/ref_tool --goldens tests/golden
      Writes the golden vectors the C test (tests/test_ribbit.c) checks: <name>.raw (reference
      TX audio or impaired channel), <name>.events (reference decoder output), index.tsv.

Needs python3 and numpy (channel impairments). Deterministic: fixed seeds.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

RATES = (8000, 16000, 32000, 44100, 48000)


def utf8_upto(text, limit):
    """The longest prefix of text whose UTF-8 encoding fits limit bytes, padded with '.'."""
    while len(text.encode()) > limit:
        text = text[:-1]
    data = text.encode()
    return data + b"." * (limit - len(data))


TEXTS = {
    "ping": b"",
    "short": b"CQ CQ DE DG6FL",
    "m16max": bytes(range(0x21, 0x21 + 85)),
    "m15min": bytes((0x41 + i % 26) for i in range(86)),
    "m15max": utf8_upto("Grüße aus Mainz, Ribbit über FM. " * 4, 128),
    "m14min": bytes((0x61 + i % 26) for i in range(129)),
    "m14max": bytes((0x20 + (i * 7) % 95) for i in range(170)),
}


def run(cmd):
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=600)
    if res.returncode:
        raise RuntimeError(f"{cmd[0]} failed: {res.stderr.strip()}")
    return res.stdout


def tx(tool, out, rate, text, call="DG6FL", carrier=1500, noise=0, fancy=0):
    payload = text.hex() if text else "-"
    run([tool, "tx", rate, carrier, noise, fancy, call, payload, out])
    return np.fromfile(out, dtype="<i2")


def shift(x, hz, rate):
    """Real signal shifted by hz (analytic signal via FFT, then mixed)."""
    n = len(x)
    spec = np.fft.fft(x.astype(np.float64))
    h = np.zeros(n)
    h[0] = 1
    h[1:(n + 1) // 2] = 2
    if n % 2 == 0:
        h[n // 2] = 1
    analytic = np.fft.ifft(spec * h)
    t = np.arange(n) / rate
    return np.real(analytic * np.exp(2j * np.pi * hz * t))


def channel(sig, rate, seed, snr_db=None, offset_hz=0.0, gain=1.0, lead_s=0.5, tail_s=0.5,
            cut=0):
    rng = np.random.default_rng(seed)
    x = np.concatenate([np.zeros(int(lead_s * rate)), sig.astype(np.float64),
                        np.zeros(int(tail_s * rate))])
    if offset_hz:
        x = shift(x, offset_hz, rate)
    x = x * gain
    if snr_db is not None:
        p = np.mean(sig.astype(np.float64) ** 2) * gain * gain
        x = x + rng.normal(0, np.sqrt(p / 10 ** (snr_db / 10)), len(x))
    if cut:
        x = x[:-cut]
    return np.clip(np.round(x), -32768, 32767).astype("<i2")


def vectors():
    """(name, rate, text, call, carrier, noise, fancy, impairments) of the differential set."""
    out = []
    for rate in RATES:
        for key, text in TEXTS.items():
            out.append((f"{rate}_{key}", rate, text, "DG6FL", 1500, 0, 0, {}))
        out.append((f"{rate}_carrier1000_noise3_fancy", rate, TEXTS["short"], "DL1ABC",
                    1000, 3, 1, {}))
        out.append((f"{rate}_carrier2000", rate, TEXTS["m15max"], "DB0XYZ", 2000, 1, 0, {}))
        for snr in (20, 10, 6, 3, 0):
            out.append((f"{rate}_snr{snr}", rate, TEXTS["m14max"], "DG6FL", 1500, 0, 0,
                        {"snr_db": snr}))
        for hz in (-50.0, -7.3, 12.5, 37.1, 50.0):
            out.append((f"{rate}_offset{hz:+g}", rate, TEXTS["short"], "DG6FL", 1500, 0, 0,
                        {"offset_hz": hz, "snr_db": 15}))
        out.append((f"{rate}_quiet", rate, TEXTS["m16max"], "DG6FL", 1500, 0, 0,
                    {"gain": 0.01}))
        out.append((f"{rate}_loud", rate, TEXTS["m16max"], "DG6FL", 1500, 0, 0,
                    {"gain": 4.0}))
        out.append((f"{rate}_truncated", rate, TEXTS["m16max"], "DG6FL", 1500, 0, 0,
                    {"tail_s": 0.0, "cut": rate // 10}))
    return out


def compare(core, ref, work):
    work.mkdir(parents=True, exist_ok=True)
    tx_ok = rx_ok = 0
    failures = []
    decoded = {}
    for i, (name, rate, text, call, carrier, noise, fancy, imp) in enumerate(vectors()):
        a = tx(ref, work / f"{name}.ref.raw", rate, text, call, carrier, noise, fancy)
        b = tx(core, work / f"{name}.core.raw", rate, text, call, carrier, noise, fancy)
        if a.tobytes() != b.tobytes():
            failures.append(f"TX {name}: {len(a)} vs {len(b)} samples, "
                            f"{int(np.sum(a[:min(len(a), len(b))] != b[:min(len(a), len(b))]))}"
                            " differ")
        else:
            tx_ok += 1
        ch = channel(a, rate, seed=1000 + i, **imp)
        path = work / f"{name}.ch.raw"
        ch.tofile(path)
        ev_ref = run([ref, "rx", rate, path])
        ev_core = run([core, "rx", rate, path])
        if ev_ref != ev_core:
            failures.append(f"RX {name}:\n  ref : {ev_ref.strip()}\n  core: {ev_core.strip()}")
        else:
            rx_ok += 1
        want = "PING" if not text else "DONE"
        decoded[name] = any(ln.startswith(want) and (not text or ln.endswith(text.hex()))
                            for ln in ev_core.splitlines())
    clean = [n for n, *_, imp in vectors() if not imp or set(imp) == {"gain"}]
    missed = [n for n in clean if not decoded[n]]
    n = len(vectors())
    print(f"TX identical: {tx_ok}/{n}")
    print(f"RX events identical: {rx_ok}/{n}")
    print("decoded: " + ", ".join(f"{k}={'y' if v else 'n'}" for k, v in decoded.items()))
    for f in failures:
        print("MISMATCH " + f)
    for m in missed:
        print(f"NOT DECODED (clean channel) {m}")
    return 1 if failures or missed else 0


GOLDEN = (
    ("tx8000_short", 8000, TEXTS["short"], "DG6FL", 1500, 0, 0, None),
    ("tx8000_ping", 8000, TEXTS["ping"], "DG6FL", 1500, 0, 0, None),
    ("tx8000_m15_fancy", 8000, TEXTS["m15max"], "DB0XYZ", 1000, 2, 1, None),
    ("tx8000_m14", 8000, TEXTS["m14max"], "DL1ABC", 1500, 0, 0, None),
    ("tx48000_m16", 48000, TEXTS["m16max"], "DG6FL", 1500, 1, 0, None),
    ("tx44100_short", 44100, TEXTS["short"], "DG6FL", 2000, 0, 0, None),
    ("rx8000_snr3", 8000, TEXTS["m14max"], "DG6FL", 1500, 0, 0, {"snr_db": 3}),
    ("rx8000_offset37", 8000, TEXTS["short"], "DG6FL", 1500, 0, 0,
     {"offset_hz": 37.1, "snr_db": 15}),
    ("rx16000_offset-7", 16000, TEXTS["m16max"], "DG6FL", 1500, 0, 0,
     {"offset_hz": -7.3, "snr_db": 10}),
    ("rx8000_noise_only", 8000, None, None, None, None, None, {"snr_db": None}),
)


def goldens(ref, out):
    out.mkdir(parents=True, exist_ok=True)
    rows = ["# name\trate\tkind\ttx args (call carrier noise fancy payload_hex)"]
    for i, (name, rate, text, call, carrier, noise, fancy, imp) in enumerate(GOLDEN):
        raw = out / f"{name}.raw"
        if text is None:
            rng = np.random.default_rng(42)
            rng.normal(0, 3000, 3 * rate).round().astype("<i2").tofile(raw)
            rows.append(f"{name}\t{rate}\trx\t-")
        else:
            sig = tx(ref, raw, rate, text, call, carrier, noise, fancy)
            if imp is None:
                rows.append(f"{name}\t{rate}\ttx\t{call} {carrier} {noise} {fancy} "
                            f"{text.hex() if text else '-'}")
                ch = channel(sig, rate, seed=7 + i)
                ch.tofile(out / f"{name}.ch.raw")
                (out / f"{name}.events").write_text(run([ref, "rx", rate, out / f"{name}.ch.raw"]))
                continue
            channel(sig, rate, seed=7 + i, **imp).tofile(raw)
            rows.append(f"{name}\t{rate}\trx\t-")
        (out / f"{name}.events").write_text(run([ref, "rx", rate, raw]))
    (out / "index.tsv").write_text("\n".join(rows) + "\n")
    print(f"wrote {len(GOLDEN)} golden vectors to {out}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--core", type=Path)
    ap.add_argument("--ref", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("build/refcheck"))
    ap.add_argument("--goldens", type=Path)
    args = ap.parse_args()
    if args.goldens:
        return goldens(args.ref, args.goldens)
    if not args.core:
        ap.error("--core is required for the differential check")
    return compare(args.core, args.ref, args.work)


if __name__ == "__main__":
    sys.exit(main())
