/*
 * ribbit-core -- Ribbit / Rattlegram COFDM text modem (receiver and transmitter) in C99.
 *
 * Copyright (c) 2026 Gerolf Ziegenhain
 *
 * A port of the Rattlegram modem by Ahmet Inan (0BSD licence,
 * https://github.com/aicodix/rattlegram, mirrored as https://github.com/wave-owl/rattlegram,
 * tag v1.19, commit 4df0f3af3c6da63f0d70512392e80c7a31b72099): app/src/main/cpp/encoder.hh,
 * decoder.hh and the aicodix DSP/code headers they include (fft.hh, schmidl_cox.hh,
 * polar*.hh, osd.hh, bose_chaudhuri_hocquenghem_encoder.hh, papr.hh, hilbert.hh, ...).
 * Ribbit is the Open Research Institute's VHF/UHF text messaging project; Rattlegram is
 * its modem. See NOTICE for the upstream notice.
 *
 * Permission to use, copy, modify, and/or distribute this software for any purpose with or
 * without fee is hereby granted.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 * THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 * SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 * ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 * CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 */
#ifndef RIBBIT_H
#define RIBBIT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RIBBIT_VERSION_STRING "0.1.0"
#define RIBBIT_ABI 1
/* Largest payload of one transmission (operation mode 14, encoder.hh: 1360 data bits). */
#define RIBBIT_MAX_PAYLOAD 170
/* Callsign: at most 9 base-37 characters (0-9, A-Z, space), encoder.hh base37(). */
#define RIBBIT_MAX_CALL 9

/* "ribbit-core 0.1.0 contract=ribbit-c1 abi=1 tx=1 list=16". Static string. The contract
 * tag changes only with an incompatible API change; list=16 is the polar list size (the
 * reference's SSE4.1/NEON build; an AVX2 build of the reference uses 32). */
const char *ribbit_version(void);
int ribbit_abi(void);                         /* RIBBIT_ABI */

#define RIBBIT_SIZEOF_EVENT 0
size_t ribbit_sizeof(int what);               /* 0 for an unknown structure */

/* The reference supports exactly these sample rates (native-lib.cpp): 8000, 16000, 32000,
 * 44100, 48000. 1 when supported. */
int ribbit_rate_supported(int rate);
/* Samples per OFDM symbol including the guard interval at this rate (1440 at 8000 Hz,
 * 160 ms); 0 for an unsupported rate. */
int ribbit_symbol_samples(int rate);

/* Operation modes (preamble meta data, decoder.hh): 0 = ping (callsign only),
 * 16 = up to 85 bytes, 15 = up to 128 bytes, 14 = up to 170 bytes. */
enum { RIBBIT_MODE_PING = 0, RIBBIT_MODE_14 = 14, RIBBIT_MODE_15 = 15, RIBBIT_MODE_16 = 16 };

/* ------------------------------------------------------------------------------------
 * Receiver: the reference Decoder<RATE> (decoder.hh) driven like Rattlegram's app
 * (MainActivity.java: feed rate/50 samples, process() whenever feed() completes a symbol,
 * fetch() on STATUS_DONE). Mono 16-bit audio, the FM-demodulated channel.
 * ------------------------------------------------------------------------------------ */
enum {
    RIBBIT_EV_SYNC = 1,   /* preamble decoded: mode 14..16, call, cfo_hz (STATUS_SYNC) */
    RIBBIT_EV_DONE = 2,   /* payload decoded, CRC good: payload, length, flips (STATUS_DONE) */
    RIBBIT_EV_FAIL = 3,   /* payload decoding failed: CRC of every list path bad (fetch -1) */
    RIBBIT_EV_PING = 4,   /* mode 0 preamble: callsign only, no payload (STATUS_PING) */
    RIBBIT_EV_NOPE = 5    /* preamble good but mode unsupported or callsign invalid */
};

typedef struct ribbit_event {
    int kind;                       /* RIBBIT_EV_* */
    int mode;                       /* operation mode from the preamble */
    float cfo_hz;                   /* carrier frequency (offset) the sync found, Hz */
    char call[RIBBIT_MAX_CALL + 1]; /* base-37 callsign, right aligned, NUL terminated */
    int flips;                      /* DONE: hard-decision bits the polar decoder corrected */
    int length;                     /* DONE: payload bytes up to the first NUL */
    unsigned char payload[RIBBIT_MAX_PAYLOAD + 1]; /* DONE: payload, NUL padded */
    long long sample;               /* input sample index at which the event was raised */
} ribbit_event;

typedef struct ribbit_rx ribbit_rx;

ribbit_rx *ribbit_rx_new(int rate);           /* NULL: unsupported rate or no memory */
void ribbit_rx_free(ribbit_rx *rx);
/* Feed n samples of any length; returns the number of events now queued. */
int ribbit_rx_process(ribbit_rx *rx, const int16_t *samples, int n);
/* Pops the oldest event into *ev; 1 when an event was returned, 0 when the queue is empty.
 * The queue holds 64 events; older events are dropped (counted in ribbit_rx_dropped). */
int ribbit_rx_next_event(ribbit_rx *rx, ribbit_event *ev);
long ribbit_rx_dropped(const ribbit_rx *rx);

/* ------------------------------------------------------------------------------------
 * Transmitter: the reference Encoder<RATE> (encoder.hh).
 * ------------------------------------------------------------------------------------ */
typedef struct ribbit_tx ribbit_tx;

ribbit_tx *ribbit_tx_new(int rate);           /* NULL: unsupported rate or no memory */
void ribbit_tx_free(ribbit_tx *tx);
/* Prepares one transmission. payload: len bytes (0 = ping), at most RIBBIT_MAX_PAYLOAD and
 * without NUL bytes; call: 1..9 characters of 0-9, A-Z, a-z or space; carrier_hz: centre of
 * the 256 payload carriers (Rattlegram's default 1500); noise_symbols: leading noise
 * symbols (0..); fancy_header: 1 = the callsign raster after the payload.
 * Returns 0, or -1 for invalid arguments. */
int ribbit_tx_configure(ribbit_tx *tx, const unsigned char *payload, int len, const char *call,
                        int carrier_hz, int noise_symbols, int fancy_header);
/* Writes ribbit_symbol_samples(rate) samples into out. Returns 1 while the transmission
 * continues, 0 once it is over (out is then silence). */
int ribbit_tx_produce(ribbit_tx *tx, int16_t *out);

#ifdef __cplusplus
}
#endif

#endif
