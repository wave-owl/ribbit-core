# ribbit-core

Ribbit / Rattlegram text modem in C99: receiver and transmitter, no dependencies besides libm.

Ribbit is the Open Research Institute's text messaging mode for VHF/UHF FM channels; its
modem is Ahmet Inan's **Rattlegram** (<https://github.com/aicodix/rattlegram>, 0BSD). One
transmission carries up to 170 bytes of UTF-8 text plus a callsign in about one second of
audio:

| part | content |
|---|---|
| optional noise symbols | 256 carriers of pseudo noise (squelch/AGC opener) |
| Schmidl-Cox symbol | 127 differentially coded BPSK carriers on even bins (timing, coarse CFO) |
| preamble symbol | 255 DBPSK carriers: BCH(255, 71) of 55 bits meta data (base-37 callsign << 8, operation mode) + CRC-16 |
| 4 payload symbols | 256 DQPSK carriers each = 2048 code bits, polar code + CRC-32 (CA-SCL, list 16) |
| optional fancy header | the callsign drawn as an 8x11 raster in the spectrum |

Symbol 1280 samples at 8 kHz (160 ms) plus a guard of 1/8; 256 carriers at 6.25 Hz spacing =
1.6 kHz occupied bandwidth around the carrier frequency (Rattlegram's default 1500 Hz).
Operation modes: 16 = up to 85 bytes (polar (2048, 712)), 15 = up to 128 bytes (2048, 1056),
14 = up to 170 bytes (2048, 1392), 0 = ping (callsign only). Sample rates 8000, 16000, 32000,
44100 and 48000 Hz, as in the reference.

## Bit-exact port

The core is a line-by-line port of the reference's `encoder.hh` and `decoder.hh` and the
aicodix headers they include (`fft.hh`, `unit_circle.hh`, `schmidl_cox.hh`, `hilbert.hh`,
`window.hh`, `kahan.hh`, `blockdc.hh`, `bip_buffer.hh`, `swa.hh`, `theil_sen.hh`,
`quick.hh`, `papr.hh`, `psk.hh`, `osd.hh`, `sort.hh`, `bose_chaudhuri_hocquenghem_encoder.hh`,
`polar*.hh`, `crc.hh`, `mls.hh`, `xorshift.hh`, `bitman.hh`). Every float expression keeps
the reference's operand order; both sides are built with `-ffp-contract=off`.

Reference: `ref/rattlegram` (submodule, mirror `wave-owl/rattlegram`, tag `v1.19`, commit
`4df0f3af3c6da63f0d70512392e80c7a31b72099`), compiled unmodified by `tools/ref_tool.cpp`.

* `make ref` renders 110 vectors (5 rates x 7 payload sizes incl. ping and every mode
  boundary, other carriers, noise symbols, fancy header, AWGN 20..0 dB, frequency offsets
  -50..+50 Hz, level x0.01 and x4, truncated) with both and decodes the channels with both:
  **TX audio byte-identical 110/110, decoder events identical 110/110** (kind, mode,
  callsign, CFO as hex float, polar bit flips, payload, sample index).
* `make test` checks 10 golden vectors written by the reference (`tests/golden`, regenerate
  with `make goldens`), round trips at every rate and payload size, independence of the
  caller's block size and the argument checks. `make sanitize` runs them under ASan/UBSan.

Measured on Apple silicon (arm64, clang 17): -O2 and -O3 builds give identical results.

### Deviations from the reference

1. The spectrum/spectrogram display (`Decoder::spectrum`) is not ported.
2. Mono input/output only (`channel_select` 0); stereo selection is the caller's business.
3. The polar list decoder has 16 paths: the reference's SSE4.1 and NEON builds. An AVX2 build
   of the reference uses 32 paths and can decode marginal frames differently.
4. Undefined behaviour of the reference is defined: a NaN soft value quantizes to 0 (the
   reference converts NaN to `int8_t`); padding of the OSD matrix and the Theil-Sen buffer
   are zero-initialised (the reference reads uninitialised memory only in a degenerate case).
5. The receiver is driven exactly like Rattlegram's app (`rate/50` samples per `feed()`,
   `process()` when `feed()` returns true, `fetch()` on DONE); the core buffers internally,
   so the caller's block size never changes the result.
6. The transmitter rejects what the app's UI never lets through: payloads over 170 bytes or
   with NUL bytes, callsigns that are empty, longer than 9 characters or contain characters
   outside 0-9, A-Z, a-z and space.
7. `cmplx(nrz(noise_seq()), nrz(noise_seq()))` evaluates the real part first, as clang (and
   the Android NDK) does; a gcc-built reference may differ in the noise symbols.
8. The reference built against libstdc++ (Linux, or gcc on macOS) reports the Schmidl-Cox
   CFO one ulp off in 2 of the 110 vectors (`44100_offset+12.5`, `48000_offset-50`; payloads,
   flips and sample indices identical). The authority is the reference built with Apple
   clang and libc++, which wrote the goldens; CI runs `make ref` there and `make test
   sanitize` (goldens) on Linux and macOS.

A complete transmission only decodes once the receiver has seen roughly two symbols past its
end (the reference processes symbols from a 4-symbol buffer); a file decode should append
about 0.5 s of silence.

## API

`include/ribbit.h`: `ribbit_version()` (`"ribbit-core 0.1.0 contract=ribbit-c1 abi=1 tx=1
list=16"`), `ribbit_abi()`, `ribbit_sizeof()`, `ribbit_rate_supported()`,
`ribbit_symbol_samples()`; receiver `ribbit_rx_new/free/process/next_event/dropped` with
`ribbit_event` (SYNC, DONE, FAIL, PING, NOPE); transmitter `ribbit_tx_new/free/configure/produce`.

## Build

```
make                 # build/libribbit.a, build/libribbit.{dylib,so}
make test sanitize   # C tests
make ref             # differential check against the C++ reference (needs ref/rattlegram, numpy)
```

## Licence

0BSD, like the reference. See LICENSE and NOTICE.

<!-- wave-owl-source-sbom-docs:start -->
## Source SBOMs

Every published release, including prereleases, receives automatically generated
CycloneDX JSON and readable documentation. See the [initial source report](SBOM.md),
[CycloneDX JSON](SBOM.cdx.json) and [release policy](SBOM-POLICY.md) for scope,
release assets, provenance declarations and reproduction instructions.
<!-- wave-owl-source-sbom-docs:end -->
