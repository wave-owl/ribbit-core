/*
 * ribbit-core internal helpers. Every function names the reference header it ports; the
 * arithmetic keeps the reference's expression order and float precision (built with
 * -ffp-contract=off on both sides), which is what makes the outputs bit-identical.
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#ifndef RIBBIT_INTERNAL_H
#define RIBBIT_INTERNAL_H

#include <math.h>
#include <stdint.h>
#include <stddef.h>

/* ---- complex.hh: DSP::Complex<float> ---- */
typedef struct { float re, im; } cf;

static inline cf cf_make(float re, float im) { cf c; c.re = re; c.im = im; return c; }
static inline cf cf_add(cf a, cf b) { return cf_make(a.re + b.re, a.im + b.im); }
static inline cf cf_sub(cf a, cf b) { return cf_make(a.re - b.re, a.im - b.im); }
static inline cf cf_scale(float a, cf b) { return cf_make(a * b.re, a * b.im); }
static inline cf cf_div_real(cf a, float b) { return cf_make(a.re / b, a.im / b); }
static inline cf cf_mul(cf a, cf b)
{
    return cf_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline cf cf_div(cf a, cf b)
{
    return cf_make((a.re * b.re + a.im * b.im) / (b.re * b.re + b.im * b.im),
                   (a.im * b.re - a.re * b.im) / (b.re * b.re + b.im * b.im));
}
static inline cf cf_conj(cf a) { return cf_make(a.re, -a.im); }
static inline float cf_norm(cf a) { return a.re * a.re + a.im * a.im; }
static inline float cf_abs(cf a) { return sqrtf(cf_norm(a)); }
static inline float cf_arg(cf a) { return atan2f(a.im, a.re); }
static inline cf cf_polar(float a, float b) { return cf_make(a * cosf(b), a * sinf(b)); }

/* ---- const.hh (the long double literals, rounded like the reference) ---- */
#define RB_PI ((float)3.141592653589793238462643383279502884197169399375105820974944592307L)
#define RB_TWO_PI ((float)6.283185307179586476925286766559005768394338798750211641949889184615L)
#define RB_INV_SQRT_TWO \
    ((float).7071067811865475244008443621048490392848359376884740365883398689953L)

/* ---- unit_circle.hh ---- */
float rb_uc_cos(int n, int N);
float rb_uc_sin(int n, int N);

/* ---- phasor.hh ---- */
typedef struct { cf prev, delta; } rb_phasor;
static inline void rb_phasor_init(rb_phasor *p) { p->prev = cf_make(1, 0); p->delta = cf_make(1, 0); }
static inline void rb_phasor_omega(rb_phasor *p, float v) { p->delta = cf_make(cosf(v), sinf(v)); }
static inline void rb_phasor_omega_n(rb_phasor *p, int n, int N)
{
    p->delta = cf_make(rb_uc_cos(n, N), rb_uc_sin(n, N));
}
static inline cf rb_phasor_next(rb_phasor *p)
{
    cf tmp = p->prev;
    p->prev = cf_mul(p->delta, p->prev);
    p->prev = cf_div_real(p->prev, cf_abs(p->prev));
    return tmp;
}

/* ---- fft.hh: mixed-radix decimation-in-time FFT, radices 2, 3, 4, 5, 7, 8 ---- */
typedef struct {
    int n, sign;
    cf *z;                       /* factors[k] = (cos(k, n), sign * sin(k, n)) */
} rb_fft;
int rb_fft_init(rb_fft *f, int n, int sign);   /* 0, or -1: size not supported */
void rb_fft_free(rb_fft *f);
void rb_fft_run(const rb_fft *f, cf *out, const cf *in);

/* ---- mls.hh ---- */
typedef struct { int poly, test, reg; } rb_mls;
void rb_mls_init(rb_mls *m, int poly);
int rb_mls_next(rb_mls *m);

/* ---- xorshift.hh: Xorshift32 ---- */
typedef struct { uint32_t y; } rb_xorshift32;
static inline void rb_xorshift32_init(rb_xorshift32 *x) { x->y = 2463534242u; }
static inline uint32_t rb_xorshift32_next(rb_xorshift32 *x)
{
    x->y ^= x->y << 13;
    x->y ^= x->y >> 17;
    x->y ^= x->y << 5;
    return x->y;
}

/* ---- bitman.hh ---- */
static inline void rb_xor_be_bit(uint8_t *buf, int pos, int val)
{
    buf[pos / 8] ^= (uint8_t)(val << (7 - pos % 8));
}
static inline void rb_set_be_bit(uint8_t *buf, int pos, int val)
{
    buf[pos / 8] = (uint8_t)((~(1 << (7 - pos % 8)) & buf[pos / 8]) | (val << (7 - pos % 8)));
}
static inline void rb_set_le_bit(uint8_t *buf, int pos, int val)
{
    buf[pos / 8] = (uint8_t)((~(1 << (pos % 8)) & buf[pos / 8]) | (val << (pos % 8)));
}
static inline int rb_get_be_bit(const uint8_t *buf, int pos) { return (buf[pos / 8] >> (7 - pos % 8)) & 1; }
static inline int rb_get_le_bit(const uint8_t *buf, int pos) { return (buf[pos / 8] >> (pos % 8)) & 1; }

/* ---- crc.hh ---- */
uint16_t rb_crc16_u64(uint64_t data);               /* CRC<uint16_t>(0xA8F4) of 8 bytes */
typedef struct { uint32_t lut[256]; uint32_t crc; } rb_crc32;
void rb_crc32_init(rb_crc32 *c);                    /* CRC<uint32_t>(0x8F6E37A0) */
static inline void rb_crc32_bit(rb_crc32 *c, int bit)
{
    uint32_t tmp = c->crc ^ (uint32_t)(bit != 0);
    c->crc = (c->crc >> 1) ^ ((tmp & 1) * 0x8F6E37A0u);
}
static inline void rb_crc32_byte(rb_crc32 *c, uint8_t data)
{
    uint32_t tmp = c->crc ^ data;
    c->crc = (c->crc >> 8) ^ c->lut[tmp & 255];
}

/* ---- bose_chaudhuri_hocquenghem_encoder.hh / osd.hh: BCH(255, 71) ---- */
void rb_bch_encode(const uint8_t *data, uint8_t *parity);       /* 71 bits -> 184 bits */
void rb_bch_matrix(int8_t *genmat);                             /* systematic 71 x 255 */
typedef struct rb_osd rb_osd;
rb_osd *rb_osd_new(void);
void rb_osd_free(rb_osd *o);
int rb_osd_decode(rb_osd *o, uint8_t *hard, const int8_t *soft, const int8_t *genmat);

/* ---- polar.hh / polar_tables.hh ---- */
extern const uint32_t rb_frozen_2048_1392[64];
extern const uint32_t rb_frozen_2048_1056[64];
extern const uint32_t rb_frozen_2048_712[64];
/* PolarEncoder<int8_t>: 2048 code bits from data_bits message bits + CRC-32. */
void rb_polar_encode(int8_t *code, const uint8_t *message, const uint32_t *frozen, int data_bits);
typedef struct rb_polar_dec rb_polar_dec;
rb_polar_dec *rb_polar_dec_new(void);
void rb_polar_dec_free(rb_polar_dec *d);
/* PolarDecoder<int8_t> with a list of 16: flips >= 0, or -1 when no path passes the CRC. */
int rb_polar_decode(rb_polar_dec *d, uint8_t *message, const int8_t *code,
                    const uint32_t *frozen, int data_bits);

/* ---- psk.hh ---- */
#define RB_RCP_SQRT_2 ((float)0.70710678118654752440)
static inline cf rb_qpsk_map(const int8_t *b)
{
    return cf_scale(RB_RCP_SQRT_2, cf_make((float)b[0], (float)b[1]));
}
static inline void rb_qpsk_hard(int8_t *b, cf c)
{
    b[0] = c.re < 0.0f ? -1 : 1;
    b[1] = c.im < 0.0f ? -1 : 1;
}
static inline int8_t rb_quantize(float dist, float precision, float value)
{
    value *= dist * precision;
    value = nearbyintf(value);
    /* std::min<float>(std::max<float>(value, -128), 127), including its NaN behaviour */
    value = value < -128.0f ? -128.0f : value;
    value = 127.0f < value ? 127.0f : value;
    return value != value ? (int8_t)0 : (int8_t)value;
}
/* PhaseShiftKeying<4>::DIST = 2 * rcp_sqrt_2 and PhaseShiftKeying<2>::DIST = 2. */
#define RB_QPSK_DIST (2 * RB_RCP_SQRT_2)
#define RB_BPSK_DIST 2.0f

/* ---- base37_bitmap.hh ---- */
extern const uint8_t rb_base37_bitmap[407];

/* Shared modem constants (encoder.hh / decoder.hh). */
enum {
    RB_CODE_ORDER = 11,
    RB_CODE_LEN = 1 << 11,
    RB_MOD_BITS = 2,
    RB_SYMBOL_COUNT = 4,
    RB_COR_SEQ_LEN = 127,
    RB_COR_SEQ_OFF = 1 - 127,
    RB_COR_SEQ_POLY = 0x89,     /* 0b10001001 */
    RB_PRE_SEQ_LEN = 255,
    RB_PRE_SEQ_OFF = -255 / 2,
    RB_PRE_SEQ_POLY = 0x12B,    /* 0b100101011 */
    RB_PAY_CAR_CNT = 256,
    RB_PAY_CAR_OFF = -256 / 2,
    RB_NOISE_POLY = 0x951       /* 0b100101010001 */
};

static inline int rb_nrz(int bit) { return 1 - 2 * (bit != 0); }
/* Frozen bit table and data bits of an operation mode; data_bits 0 for ping/unknown. */
const uint32_t *rb_mode_frozen(int mode, int *data_bits);

#endif
