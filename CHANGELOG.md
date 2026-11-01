# Changelog

## Unreleased

- CI (GitHub Actions, Linux + macOS): `make`, `make test`, `make sanitize` (ASan/UBSan) and
  the differential check `make ref` against the C++ reference (macOS, Apple clang + libc++).
- `tools/ref_tool.cpp` includes `<cstdint>` first, so the reference builds with libstdc++.
- README: a libstdc++ build of the reference reports the CFO one ulp off in 2 of 110 vectors.

## 0.1.0

- Receiver and transmitter ported from Rattlegram v1.19 (`decoder.hh`, `encoder.hh`), contract
  `ribbit-c1`, ABI 1, `tx=1`, polar list 16.
- Differential check against the unmodified C++ reference (`make ref`): 110/110 TX vectors
  byte-identical, 110/110 decoder event lists identical; 10 golden vectors in `tests/golden`.
