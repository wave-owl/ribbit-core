# Changelog

## 0.1.0

- Receiver and transmitter ported from Rattlegram v1.19 (`decoder.hh`, `encoder.hh`), contract
  `ribbit-c1`, ABI 1, `tx=1`, polar list 16.
- Differential check against the unmodified C++ reference (`make ref`): 110/110 TX vectors
  byte-identical, 110/110 decoder event lists identical; 10 golden vectors in `tests/golden`.
