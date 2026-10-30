/*
 * Transmitter: port of Rattlegram's Encoder<RATE> (encoder.hh) with papr.hh
 * (ImprovePAPR), polar.hh (PolarEncoder) and the BCH(255, 71) preamble code. Mono output
 * (channel_select 0 of the reference).
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <stdlib.h>
#include <string.h>
#include "ribbit.h"
#include "ribbit_internal.h"

enum { MAX_BITS = 1360 };

struct ribbit_tx {
    int rate, S, G, E, fact;
    rb_fft bwd, papr_fwd, papr_bwd;
    cf *temp, *freq, prev[RB_PAY_CAR_CNT], *guard;
    cf *papr_temp, *papr_over;
    unsigned char *used;
    uint8_t mesg[MAX_BITS / 8], call[9];
    int8_t code[RB_CODE_LEN];
    uint64_t meta_data;
    rb_mls noise_seq;
    int operation_mode, carrier_offset, symbol_number, count_down, fancy_line, noise_count;
};

static int bin(const ribbit_tx *tx, int carrier)
{
    return (carrier + tx->carrier_offset + tx->S) % tx->S;
}

/* papr.hh: ImprovePAPR<cmplx, size, fact> */
static void improve_papr(ribbit_tx *tx, cf *freq)
{
    const int size = tx->S, fact = tx->fact, n = fact * size;
    float factor;
    int i;
    for (i = 0; i < size; ++i)
        tx->used[i] = freq[i].re != 0 || freq[i].im != 0;
    if (fact == 1) {
        rb_fft_run(&tx->papr_bwd, tx->papr_temp, freq);
    } else {
        for (i = 0; i < size / 2; ++i)
            tx->papr_over[i] = freq[i];
        for (i = size / 2; i < n - size / 2; ++i)
            tx->papr_over[i] = cf_make(0, 0);
        for (i = size / 2; i < size; ++i)
            tx->papr_over[size * (fact - 1) + i] = freq[i];
        rb_fft_run(&tx->papr_bwd, tx->papr_temp, tx->papr_over);
    }
    factor = 1 / sqrtf((float)n);
    for (i = 0; i < n; ++i)
        tx->papr_temp[i] = cf_scale(factor, tx->papr_temp[i]);
    for (i = 0; i < n; ++i) {
        float pwr = cf_norm(tx->papr_temp[i]);
        if (pwr > 1)
            tx->papr_temp[i] = cf_div_real(tx->papr_temp[i], sqrtf(pwr));
    }
    if (fact == 1) {
        rb_fft_run(&tx->papr_fwd, freq, tx->papr_temp);
        for (i = 0; i < size; ++i)
            freq[i] = tx->used[i] ? cf_scale(factor, freq[i]) : cf_make(0, 0);
    } else {
        rb_fft_run(&tx->papr_fwd, tx->papr_over, tx->papr_temp);
        for (i = 0; i < size / 2; ++i)
            if (tx->used[i])
                freq[i] = cf_scale(factor, tx->papr_over[i]);
        for (i = size / 2; i < size; ++i)
            if (tx->used[i])
                freq[i] = cf_scale(factor, tx->papr_over[size * (fact - 1) + i]);
    }
}

static void transform(ribbit_tx *tx)
{
    int i;
    improve_papr(tx, tx->freq);
    rb_fft_run(&tx->bwd, tx->temp, tx->freq);
    for (i = 0; i < tx->S; ++i)
        tx->temp[i] = cf_div_real(tx->temp[i], sqrtf((float)(8 * tx->S)));
}

static void clear_freq(ribbit_tx *tx)
{
    int i;
    for (i = 0; i < tx->S; ++i)
        tx->freq[i] = cf_make(0, 0);
}

static void schmidl_cox(ribbit_tx *tx)
{
    rb_mls seq;
    float factor = sqrtf((float)(2 * tx->S) / RB_COR_SEQ_LEN);
    int i;
    rb_mls_init(&seq, RB_COR_SEQ_POLY);
    clear_freq(tx);
    tx->freq[bin(tx, RB_COR_SEQ_OFF - 2)] = cf_make(factor, 0);
    for (i = 0; i < RB_COR_SEQ_LEN; ++i)
        tx->freq[bin(tx, 2 * i + RB_COR_SEQ_OFF)] = cf_make((float)rb_nrz(rb_mls_next(&seq)), 0);
    for (i = 0; i < RB_COR_SEQ_LEN; ++i)
        tx->freq[bin(tx, 2 * i + RB_COR_SEQ_OFF)] =
            cf_mul(tx->freq[bin(tx, 2 * (i - 1) + RB_COR_SEQ_OFF)],
                   tx->freq[bin(tx, 2 * i + RB_COR_SEQ_OFF)]);
    transform(tx);
}

static void preamble(ribbit_tx *tx)
{
    uint8_t data[9], parity[23];
    rb_mls seq;
    float factor;
    uint16_t cs;
    int i;
    memset(data, 0, sizeof(data));
    memset(parity, 0, sizeof(parity));
    for (i = 0; i < 55; ++i)
        rb_set_be_bit(data, i, (int)((tx->meta_data >> i) & 1));
    cs = rb_crc16_u64(tx->meta_data << 9);
    for (i = 0; i < 16; ++i)
        rb_set_be_bit(data, i + 55, (cs >> i) & 1);
    rb_bch_encode(data, parity);
    rb_mls_init(&seq, RB_PRE_SEQ_POLY);
    factor = sqrtf((float)tx->S / RB_PRE_SEQ_LEN);
    clear_freq(tx);
    tx->freq[bin(tx, RB_PRE_SEQ_OFF - 1)] = cf_make(factor, 0);
    for (i = 0; i < 71; ++i)
        tx->freq[bin(tx, i + RB_PRE_SEQ_OFF)] = cf_make((float)rb_nrz(rb_get_be_bit(data, i)), 0);
    for (i = 71; i < RB_PRE_SEQ_LEN; ++i)
        tx->freq[bin(tx, i + RB_PRE_SEQ_OFF)] =
            cf_make((float)rb_nrz(rb_get_be_bit(parity, i - 71)), 0);
    for (i = 0; i < RB_PRE_SEQ_LEN; ++i)
        tx->freq[bin(tx, i + RB_PRE_SEQ_OFF)] =
            cf_mul(tx->freq[bin(tx, i - 1 + RB_PRE_SEQ_OFF)], tx->freq[bin(tx, i + RB_PRE_SEQ_OFF)]);
    for (i = 0; i < RB_PRE_SEQ_LEN; ++i)
        tx->freq[bin(tx, i + RB_PRE_SEQ_OFF)] =
            cf_scale((float)rb_nrz(rb_mls_next(&seq)), tx->freq[bin(tx, i + RB_PRE_SEQ_OFF)]);
    for (i = 0; i < RB_PAY_CAR_CNT; ++i)
        tx->prev[i] = tx->freq[bin(tx, i + RB_PAY_CAR_OFF)];
    transform(tx);
}

static void fancy_symbol(ribbit_tx *tx)
{
    const int fancy_off = -(8 * 9 * 3) / 2;
    int active_carriers = 1, i, j;
    float factor;
    for (j = 0; j < 9; ++j)
        for (i = 0; i < 8; ++i)
            active_carriers += (rb_base37_bitmap[tx->call[j] + 37 * tx->fancy_line] >> i) & 1;
    factor = sqrtf((float)tx->S / (float)active_carriers);
    clear_freq(tx);
    for (j = 0; j < 9; ++j)
        for (i = 0; i < 8; ++i)
            if (rb_base37_bitmap[tx->call[j] + 37 * tx->fancy_line] & (1 << (7 - i)))
                tx->freq[bin(tx, (8 * j + i) * 3 + fancy_off)] =
                    cf_make(factor * (float)rb_nrz(rb_mls_next(&tx->noise_seq)), 0);
    transform(tx);
}

static void noise_symbol(ribbit_tx *tx)
{
    float factor = sqrtf((float)tx->S / (float)RB_PAY_CAR_CNT);
    int i;
    clear_freq(tx);
    for (i = 0; i < RB_PAY_CAR_CNT; ++i) {
        /* cmplx(nrz(noise_seq()), nrz(noise_seq())): arguments evaluated left to right, as
         * the clang-built reference does */
        float re = (float)rb_nrz(rb_mls_next(&tx->noise_seq));
        float im = (float)rb_nrz(rb_mls_next(&tx->noise_seq));
        tx->freq[bin(tx, i + RB_PAY_CAR_OFF)] = cf_scale(factor, cf_make(re, im));
    }
    transform(tx);
}

static void payload_symbol(ribbit_tx *tx)
{
    int i;
    clear_freq(tx);
    for (i = 0; i < RB_PAY_CAR_CNT; ++i) {
        tx->prev[i] = cf_mul(rb_qpsk_map(tx->code + RB_MOD_BITS *
                                         (RB_PAY_CAR_CNT * tx->symbol_number + i)),
                             tx->prev[i]);
        tx->freq[bin(tx, i + RB_PAY_CAR_OFF)] = tx->prev[i];
    }
    transform(tx);
}

static int16_t sample16(float v)
{
    float x = nearbyintf(32767 * v);
    x = x < -32768.0f ? -32768.0f : x;   /* std::clamp<float>(x, -32768, 32767) */
    x = 32767.0f < x ? 32767.0f : x;
    return (int16_t)x;
}

int ribbit_tx_produce(ribbit_tx *tx, int16_t *out)
{
    int data_symbol = 0, i;
    if (!tx || !out)
        return 0;
    switch (tx->count_down) {
    case 5:
        if (tx->noise_count) {
            --tx->noise_count;
            noise_symbol(tx);
            break;
        }
        --tx->count_down;
        /* fall through */
    case 4:
        schmidl_cox(tx);
        data_symbol = 1;
        --tx->count_down;
        break;
    case 3:
        preamble(tx);
        data_symbol = 1;
        --tx->count_down;
        if (!tx->operation_mode)
            --tx->count_down;
        break;
    case 2:
        payload_symbol(tx);
        data_symbol = 1;
        if (++tx->symbol_number == RB_SYMBOL_COUNT)
            --tx->count_down;
        break;
    case 1:
        if (tx->fancy_line) {
            --tx->fancy_line;
            fancy_symbol(tx);
            break;
        }
        for (i = 0; i < tx->S; ++i)
            tx->temp[i] = cf_make(0, 0);
        --tx->count_down;
        break;
    default:
        for (i = 0; i < tx->E; ++i)
            out[i] = 0;
        return 0;
    }
    for (i = 0; i < tx->G; ++i) {
        float x = (float)i / (float)(tx->G - 1);
        const float ratio = 0.5f;
        float y;
        cf sum;
        if (data_symbol)
            x = (ratio < x ? ratio : x) / ratio;
        y = 0.5f * (1 - cosf(RB_PI * x));
        sum = cf_add(cf_scale(1.0f - y, tx->guard[i]), cf_scale(y, tx->temp[i + tx->S - tx->G]));
        out[i] = sample16(sum.re);
    }
    for (i = 0; i < tx->G; ++i)
        tx->guard[i] = tx->temp[i];
    for (i = 0; i < tx->S; ++i)
        out[i + tx->G] = sample16(tx->temp[i].re);
    return 1;
}

static uint8_t base37_map(char c)
{
    if (c >= '0' && c <= '9')
        return (uint8_t)(c - '0' + 1);
    if (c >= 'a' && c <= 'z')
        return (uint8_t)(c - 'a' + 11);
    if (c >= 'A' && c <= 'Z')
        return (uint8_t)(c - 'A' + 11);
    return 0;
}

int ribbit_tx_configure(ribbit_tx *tx, const unsigned char *payload, int len, const char *call,
                        int carrier_hz, int noise_symbols, int fancy_header)
{
    unsigned char buf[RIBBIT_MAX_PAYLOAD + 1];
    uint64_t acc = 0;
    const uint32_t *frozen_bits;
    int data_bits, i, n, clen;
    rb_xorshift32 scrambler;
    if (!tx || !call || len < 0 || len > RIBBIT_MAX_PAYLOAD || (len && !payload) ||
        noise_symbols < 0)
        return -1;
    clen = (int)strlen(call);
    if (clen < 1 || clen > RIBBIT_MAX_CALL)
        return -1;
    for (i = 0; i < clen; ++i)
        if (call[i] != ' ' && !base37_map(call[i]))
            return -1;
    memset(buf, 0, sizeof(buf));
    for (i = 0; i < len; ++i) {
        if (!payload[i])
            return -1;
        buf[i] = payload[i];
    }
    /* Encoder::configure */
    n = 0;
    while (n <= 128 && buf[n])
        ++n;
    if (!n)
        tx->operation_mode = 0;
    else if (n <= 85)
        tx->operation_mode = 16;
    else if (n <= 128)
        tx->operation_mode = 15;
    else
        tx->operation_mode = 14;
    tx->carrier_offset = (carrier_hz * tx->S) / tx->rate;
    for (i = 0; i < clen; ++i)
        acc = 37 * acc + base37_map(call[i]);
    tx->meta_data = (acc << 8) | (uint64_t)tx->operation_mode;
    for (i = 0; i < 9; ++i)
        tx->call[i] = 0;
    for (i = 0; i < 9 && i < clen; ++i)
        tx->call[i] = base37_map(call[i]);
    tx->symbol_number = 0;
    tx->count_down = 5;
    tx->fancy_line = 11 * (fancy_header != 0);
    tx->noise_count = noise_symbols;
    for (i = 0; i < tx->G; ++i)
        tx->guard[i] = cf_make(0, 0);
    frozen_bits = rb_mode_frozen(tx->operation_mode, &data_bits);
    if (!frozen_bits)
        return 0;
    rb_xorshift32_init(&scrambler);
    for (i = 0; i < data_bits / 8; ++i)
        tx->mesg[i] = (uint8_t)(buf[i] ^ rb_xorshift32_next(&scrambler));
    rb_polar_encode(tx->code, tx->mesg, frozen_bits, data_bits);
    return 0;
}

ribbit_tx *ribbit_tx_new(int rate)
{
    ribbit_tx *tx;
    if (!ribbit_rate_supported(rate))
        return NULL;
    tx = (ribbit_tx *)calloc(1, sizeof(ribbit_tx));
    if (!tx)
        return NULL;
    tx->rate = rate;
    tx->S = (1280 * rate) / 8000;
    tx->G = tx->S / 8;
    tx->E = tx->S + tx->G;
    tx->fact = (32000 + rate / 2) / rate;
    tx->temp = (cf *)calloc((size_t)tx->E, sizeof(cf));
    tx->freq = (cf *)calloc((size_t)tx->S, sizeof(cf));
    tx->guard = (cf *)calloc((size_t)tx->G, sizeof(cf));
    tx->papr_temp = (cf *)calloc((size_t)(tx->fact * tx->S), sizeof(cf));
    tx->papr_over = (cf *)calloc((size_t)(tx->fact * tx->S), sizeof(cf));
    tx->used = (unsigned char *)calloc((size_t)tx->S, 1);
    if (!tx->temp || !tx->freq || !tx->guard || !tx->papr_temp || !tx->papr_over || !tx->used ||
        rb_fft_init(&tx->bwd, tx->S, 1) || rb_fft_init(&tx->papr_fwd, tx->fact * tx->S, -1) ||
        rb_fft_init(&tx->papr_bwd, tx->fact * tx->S, 1)) {
        ribbit_tx_free(tx);
        return NULL;
    }
    rb_mls_init(&tx->noise_seq, RB_NOISE_POLY);
    tx->symbol_number = RB_SYMBOL_COUNT;
    return tx;
}

void ribbit_tx_free(ribbit_tx *tx)
{
    if (!tx)
        return;
    rb_fft_free(&tx->bwd);
    rb_fft_free(&tx->papr_fwd);
    rb_fft_free(&tx->papr_bwd);
    free(tx->temp);
    free(tx->freq);
    free(tx->guard);
    free(tx->papr_temp);
    free(tx->papr_over);
    free(tx->used);
    free(tx);
}
