/*
 * Receiver: port of Rattlegram's Decoder<RATE> (decoder.hh) with schmidl_cox.hh, hilbert.hh
 * (Kaiser window from window.hh, kahan.hh), blockdc.hh, bip_buffer.hh, sma.hh/swa.hh,
 * delay.hh, trigger.hh, theil_sen.hh/quick.hh and phasor.hh. The spectrum/spectrogram
 * display of the reference is left out. ribbit_rx_process() drives the decoder like
 * Rattlegram's app (MainActivity.java): blocks of rate/50 samples go to feed(); whenever
 * feed() completes a symbol, process() runs, and STATUS_DONE is followed by fetch().
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <stdlib.h>
#include <string.h>
#include "ribbit.h"
#include "ribbit_internal.h"

enum { ST_OKAY = 0, ST_FAIL = 1, ST_SYNC = 2, ST_DONE = 3, ST_NOPE = 5, ST_PING = 6 };
enum { QUEUE = 64 };

/* swa.hh: sliding window accelerator with std::plus (SMA4<..., NORM=false>) */
typedef struct { float *tree; int num, leaf; } swa_f;
typedef struct { cf *tree; int num, leaf; } swa_c;

static int swa_f_init(swa_f *s, int num)
{
    s->num = num;
    s->leaf = num;
    s->tree = (float *)calloc((size_t)(2 * num), sizeof(float));
    return s->tree ? 0 : -1;
}

static float swa_f_push(swa_f *s, float input)
{
    int child, parent;
    s->tree[s->leaf] = input;
    for (child = s->leaf, parent = s->leaf / 2; parent; child = parent, parent /= 2)
        s->tree[parent] = s->tree[child] + s->tree[child ^ 1];
    if (++s->leaf >= 2 * s->num)
        s->leaf = s->num;
    return s->tree[1];
}

static int swa_c_init(swa_c *s, int num)
{
    s->num = num;
    s->leaf = num;
    s->tree = (cf *)calloc((size_t)(2 * num), sizeof(cf));
    return s->tree ? 0 : -1;
}

static cf swa_c_push(swa_c *s, cf input)
{
    int child, parent;
    s->tree[s->leaf] = input;
    for (child = s->leaf, parent = s->leaf / 2; parent; child = parent, parent /= 2)
        s->tree[parent] = cf_add(s->tree[child], s->tree[child ^ 1]);
    if (++s->leaf >= 2 * s->num)
        s->leaf = s->num;
    return s->tree[1];
}

/* schmidl_cox.hh: SchmidlCox<float, cmplx, search_pos, symbol_len, guard_len> */
typedef struct {
    int search_pos, symbol_len, guard_len, match_len, match_del;
    rb_fft fwd, bwd;
    swa_c cor;
    swa_f pwr, match;
    float *align;               /* Delay<value, match_del> */
    int align_pos;
    float thr_low, thr_high;    /* SchmittTrigger */
    int thr_prev, fall_prev;    /* SchmittTrigger, FallingEdgeTrigger */
    cf *tmp0, *tmp1, *kern;
    float timing_max, phase_max;
    int index_max;
    int symbol_pos;
    float cfo_rad, frac_cfo;
} schmidl_cox;

static int sc_bin(const schmidl_cox *sc, int carrier)
{
    return (carrier + sc->symbol_len) % sc->symbol_len;
}

static cf sc_demod_or_erase(cf curr, cf prev, float pwr)
{
    if (cf_norm(curr) > pwr && cf_norm(prev) > pwr) {
        cf cons = cf_div(curr, prev);
        if (cf_norm(cons) < 4.0f)
            return cons;
    }
    return cf_make(0, 0);
}

static int sc_init(schmidl_cox *sc, int search_pos, int symbol_len, int guard_len, const cf *seq)
{
    int i;
    sc->search_pos = search_pos;
    sc->symbol_len = symbol_len;
    sc->guard_len = guard_len;
    sc->match_len = guard_len | 1;
    sc->match_del = (sc->match_len - 1) / 2;
    if (rb_fft_init(&sc->fwd, symbol_len, -1) || rb_fft_init(&sc->bwd, symbol_len, 1))
        return -1;
    if (swa_c_init(&sc->cor, symbol_len) || swa_f_init(&sc->pwr, 2 * symbol_len) ||
        swa_f_init(&sc->match, sc->match_len))
        return -1;
    sc->align = (float *)calloc((size_t)(sc->match_del > 0 ? sc->match_del : 1), sizeof(float));
    sc->tmp0 = (cf *)calloc((size_t)symbol_len, sizeof(cf));
    sc->tmp1 = (cf *)calloc((size_t)symbol_len, sizeof(cf));
    sc->kern = (cf *)calloc((size_t)symbol_len, sizeof(cf));
    if (!sc->align || !sc->tmp0 || !sc->tmp1 || !sc->kern)
        return -1;
    sc->align_pos = 0;
    sc->thr_low = (float)(0.17 * sc->match_len);
    sc->thr_high = (float)(0.19 * sc->match_len);
    sc->thr_prev = 0;
    sc->fall_prev = 0;
    sc->timing_max = 0;
    sc->phase_max = 0;
    sc->index_max = 0;
    sc->symbol_pos = 0;
    sc->cfo_rad = 0;
    sc->frac_cfo = 0;
    rb_fft_run(&sc->fwd, sc->kern, seq);
    for (i = 0; i < symbol_len; ++i)
        sc->kern[i] = cf_div_real(cf_conj(sc->kern[i]), (float)symbol_len);
    return 0;
}

static void sc_free(schmidl_cox *sc)
{
    rb_fft_free(&sc->fwd);
    rb_fft_free(&sc->bwd);
    free(sc->cor.tree);
    free(sc->pwr.tree);
    free(sc->match.tree);
    free(sc->align);
    free(sc->tmp0);
    free(sc->tmp1);
    free(sc->kern);
}

static int sc_run(schmidl_cox *sc, const cf *samples)
{
    const int symbol_len = sc->symbol_len;
    cf P = swa_c_push(&sc->cor, cf_mul(samples[sc->search_pos + symbol_len],
                                       cf_conj(samples[sc->search_pos + 2 * symbol_len])));
    float R = 0.5f * swa_f_push(&sc->pwr, cf_norm(samples[sc->search_pos + 2 * symbol_len]));
    float min_R = (float)(0.00001 * symbol_len);
    float timing, phase, min_pwr, peak, next;
    int collect, process, i, shift, pos_err;
    rb_phasor osc;
    R = R < min_R ? min_R : R;     /* std::max(R, min_R) */
    timing = swa_f_push(&sc->match, cf_norm(P) / (R * R));
    if (sc->match_del > 0) {
        float tmp = sc->align[sc->align_pos];
        sc->align[sc->align_pos] = cf_arg(P);
        if (++sc->align_pos >= sc->match_del)
            sc->align_pos = 0;
        phase = tmp;
    } else {
        phase = cf_arg(P);
    }

    if (sc->thr_prev) {
        if (timing < sc->thr_low)
            sc->thr_prev = 0;
    } else {
        if (timing > sc->thr_high)
            sc->thr_prev = 1;
    }
    collect = sc->thr_prev;
    process = sc->fall_prev && !collect;
    sc->fall_prev = collect;

    if (!collect && !process)
        return 0;

    if (sc->timing_max < timing) {
        sc->timing_max = timing;
        sc->phase_max = phase;
        sc->index_max = sc->match_del;
    } else if (sc->index_max < symbol_len + sc->guard_len + sc->match_del) {
        ++sc->index_max;
    }

    if (!process)
        return 0;

    sc->frac_cfo = sc->phase_max / (float)symbol_len;

    rb_phasor_init(&osc);
    rb_phasor_omega(&osc, sc->frac_cfo);
    sc->symbol_pos = sc->search_pos - sc->index_max;
    sc->index_max = 0;
    sc->timing_max = 0;
    for (i = 0; i < symbol_len; ++i)
        sc->tmp1[i] = cf_mul(samples[i + sc->symbol_pos + symbol_len], rb_phasor_next(&osc));
    rb_fft_run(&sc->fwd, sc->tmp0, sc->tmp1);
    min_pwr = 0;
    for (i = 0; i < symbol_len; ++i)
        min_pwr += cf_norm(sc->tmp0[i]);
    min_pwr /= (float)symbol_len;
    for (i = 0; i < symbol_len; ++i)
        sc->tmp1[i] = sc_demod_or_erase(sc->tmp0[i], sc->tmp0[sc_bin(sc, i - 1)], min_pwr);
    rb_fft_run(&sc->fwd, sc->tmp0, sc->tmp1);
    for (i = 0; i < symbol_len; ++i)
        sc->tmp0[i] = cf_mul(sc->kern[i], sc->tmp0[i]);
    rb_fft_run(&sc->bwd, sc->tmp1, sc->tmp0);

    shift = 0;
    peak = 0;
    next = 0;
    for (i = 0; i < symbol_len; ++i) {
        float power = cf_norm(sc->tmp1[i]);
        if (power > peak) {
            next = peak;
            peak = power;
            shift = i;
        } else if (power > next) {
            next = power;
        }
    }
    if (peak <= next * 4)
        return 0;

    pos_err = (int)nearbyintf(cf_arg(sc->tmp1[shift]) * (float)symbol_len / RB_TWO_PI);
    if (abs(pos_err) > sc->guard_len / 2)
        return 0;
    sc->symbol_pos -= pos_err;

    sc->cfo_rad = (float)shift * (RB_TWO_PI / (float)symbol_len) - sc->frac_cfo;
    if (sc->cfo_rad >= RB_PI)
        sc->cfo_rad -= RB_TWO_PI;
    return 1;
}

/* window.hh Kaiser with kahan.hh's i0 */
static float kaiser_i0(float x)
{
    float high = 1.0f, low = 0.0f, val = 1.0f;
    int n;
    for (n = 1; n < 35; ++n) {
        float in, tmp, sum, nlow;
        val *= x / (float)(2 * n);
        in = val * val;
        tmp = in - low;
        sum = high + tmp;
        nlow = (sum - high) - tmp;
        if (nlow == low && sum == high)
            return sum;
        low = nlow;
        high = sum;
    }
    return high;
}

static float kaiser(float a, int n, int N)
{
    float t = (float)(2 * n) / (float)(N - 1) - 1.0f;
    return kaiser_i0(RB_PI * a * sqrtf(1.0f - t * t)) / kaiser_i0(RB_PI * a);
}

struct ribbit_rx {
    int rate, S, G, E, taps, chunk, buffer_length, search_position;
    /* hilbert.hh */
    float *hreal, *imco, reco;
    /* blockdc.hh */
    float dc_x1, dc_y1, dc_a, dc_b;
    /* bip_buffer.hh */
    cf *bip;
    int pos0, pos1;
    int buf_off;                 /* the reference's `buf` pointer, as an offset into bip */
    rb_fft fwd;
    schmidl_cox correlator;
    rb_osd *osd;
    rb_polar_dec *polar;
    int8_t *generator;
    rb_phasor osc;
    cf *temp, *freq, prev[RB_PAY_CAR_CNT], cons[RB_PAY_CAR_CNT];
    float index[RB_PAY_CAR_CNT], phase[RB_PAY_CAR_CNT];
    float *tse_temp;             /* TheilSenEstimator temp_ */
    float tse_slope, tse_yint;
    int8_t code[RB_CODE_LEN];
    int8_t soft[RB_PRE_SEQ_LEN];
    uint8_t data[(RB_PRE_SEQ_LEN + 7) / 8];
    int symbol_number, symbol_position, stored_position, staged_position;
    int staged_mode, operation_mode, accumulated;
    float stored_cfo_rad, staged_cfo_rad;
    uint64_t staged_call;
    int stored_check, staged_check;
    /* driver state */
    int16_t *pending;
    int npending;
    long long consumed;
    float sync_cfo_hz;
    char sync_call[RIBBIT_MAX_CALL + 1];
    ribbit_event queue[QUEUE];
    int qhead, qcount;
    long dropped;
};

static cf analytic(ribbit_rx *rx, float real)
{
    /* blockdc.hh then hilbert.hh */
    float y0 = rx->dc_b * (real - rx->dc_x1) + rx->dc_a * rx->dc_y1;
    int T = rx->taps, c = (T - 1) / 2, i;
    float re, im;
    rx->dc_x1 = real;
    rx->dc_y1 = y0;
    re = rx->reco * rx->hreal[c];
    im = rx->imco[0] * (rx->hreal[c - 1] - rx->hreal[c + 1]);
    for (i = 1; i < (T - 1) / 4; ++i)
        im += rx->imco[i] * (rx->hreal[c - (2 * i + 1)] - rx->hreal[c + (2 * i + 1)]);
    memmove(rx->hreal, rx->hreal + 1, sizeof(float) * (size_t)(T - 1));
    rx->hreal[T - 1] = y0;
    return cf_make(re, im);
}

static const cf *bip_push(ribbit_rx *rx, cf input)
{
    int num = rx->buffer_length;
    rx->bip[rx->pos0] = rx->bip[rx->pos1] = input;
    if (++rx->pos0 >= 2 * num)
        rx->pos0 = 0;
    if (++rx->pos1 >= 2 * num)
        rx->pos1 = 0;
    return rx->bip + (rx->pos0 < rx->pos1 ? rx->pos0 : rx->pos1);
}

static int bin(const ribbit_rx *rx, int carrier) { return (carrier + rx->S) % rx->S; }

static cf demod_or_erase(cf curr, cf prev)
{
    cf cons;
    if (cf_norm(prev) <= 0)
        return cf_make(0, 0);
    cons = cf_div(curr, prev);
    if (cf_norm(cons) > 4)
        return cf_make(0, 0);
    return cons;
}

/* quick.hh */
static void q_swap(float *a, int i, int j) { float t = a[i]; a[i] = a[j]; a[j] = t; }

static void q_median(float *A, int a, int b, int c, int d, int e)
{
    if (A[c] < A[a])
        q_swap(A, a, c);
    if (A[d] < A[b])
        q_swap(A, b, d);
    if (A[d] < A[c]) {
        q_swap(A, c, d);
        q_swap(A, a, b);
    }
    if (A[e] < A[b])
        q_swap(A, b, e);
    if (A[e] < A[c]) {
        q_swap(A, c, e);
        if (A[c] < A[a])
            q_swap(A, a, c);
    } else if (A[c] < A[b]) {
        q_swap(A, b, c);
    }
}

static void q_partition(float *a, int *l, int *h)
{
    int half = (*h - *l) / 2, quarter = (*h - *l) / 4, middle = *l + half, i;
    float pivot;
    q_median(a, *l, *l + quarter, middle, middle + quarter, *h);
    pivot = a[middle];
    for (i = *l; i <= *h;) {
        if (a[i] < pivot)
            q_swap(a, i++, (*l)++);
        else if (a[i] > pivot)
            q_swap(a, i, (*h)--);
        else
            ++i;
    }
}

static void q_insertion(float *a, int l, int h)
{
    int i, j;
    for (i = l + 1; i <= h; ++i) {
        float t = a[i];
        for (j = i; j > l && t < a[j - 1]; --j)
            a[j] = a[j - 1];
        a[j] = t;
    }
}

static float quick_select(float *a, int k, int n)
{
    int l = 0, h = n - 1;
    while (l < h) {
        if (h - l < 32) {
            q_insertion(a, l, h);
            break;
        } else {
            int lt = l, gt = h;
            q_partition(a, &lt, &gt);
            if (k < lt)
                h = lt - 1;
            else if (k > gt)
                l = gt + 1;
            else
                break;
        }
    }
    return a[k];
}

static void tse_compute(ribbit_rx *rx, const float *x, const float *y, int len)
{
    const int size = ((RB_PAY_CAR_CNT - 1) * RB_PAY_CAR_CNT) / 2;
    int count = 0, i, j;
    for (i = 0; count < size && i < len; ++i)
        for (j = i + 1; count < size && j < len; ++j)
            if (x[j] != x[i])
                rx->tse_temp[count++] = (y[j] - y[i]) / (x[j] - x[i]);
    rx->tse_slope = quick_select(rx->tse_temp, count / 2, count);
    count = 0;
    for (i = 0; count < size && i < len; ++i)
        rx->tse_temp[count++] = y[i] - rx->tse_slope * x[i];
    rx->tse_yint = quick_select(rx->tse_temp, count / 2, count);
}

static void compensate(ribbit_rx *rx)
{
    int count = 0, i;
    for (i = 0; i < RB_PAY_CAR_CNT; ++i) {
        cf con = rx->cons[i];
        if (con.re != 0 && con.im != 0) {
            int8_t tmp[RB_MOD_BITS];
            rb_qpsk_hard(tmp, con);
            rx->index[count] = (float)(i + RB_PAY_CAR_OFF);
            rx->phase[count] = cf_arg(cf_mul(con, cf_conj(rb_qpsk_map(tmp))));
            ++count;
        }
    }
    tse_compute(rx, rx->index, rx->phase, count);
    for (i = 0; i < RB_PAY_CAR_CNT; ++i)
        rx->cons[i] = cf_mul(rx->cons[i],
                             cf_polar(1.0f, -(rx->tse_yint + rx->tse_slope *
                                              (float)(i + RB_PAY_CAR_OFF))));
}

static float precision(const ribbit_rx *rx)
{
    float sp = 0, np = 0;
    int i;
    for (i = 0; i < RB_PAY_CAR_CNT; ++i) {
        int8_t tmp[RB_MOD_BITS];
        cf hard, error;
        rb_qpsk_hard(tmp, rx->cons[i]);
        hard = rb_qpsk_map(tmp);
        error = cf_sub(rx->cons[i], hard);
        sp += cf_norm(hard);
        np += cf_norm(error);
    }
    return sp / np;
}

static void demap(ribbit_rx *rx)
{
    float pre = precision(rx);
    int i;
    for (i = 0; i < RB_PAY_CAR_CNT; ++i) {
        int8_t *b = rx->code + RB_MOD_BITS * (rx->symbol_number * RB_PAY_CAR_CNT + i);
        b[0] = rb_quantize(RB_QPSK_DIST, pre, rx->cons[i].re);
        b[1] = rb_quantize(RB_QPSK_DIST, pre, rx->cons[i].im);
    }
}

static int preamble(ribbit_rx *rx)
{
    rb_phasor nco;
    rb_mls seq;
    uint64_t md = 0;
    uint16_t cs = 0;
    int i;
    const cf *buf = rx->bip + rx->buf_off;
    rb_phasor_init(&nco);
    rb_phasor_omega(&nco, -rx->staged_cfo_rad);
    for (i = 0; i < rx->S; ++i)
        rx->temp[i] = cf_mul(buf[rx->staged_position + i], rb_phasor_next(&nco));
    rb_fft_run(&rx->fwd, rx->freq, rx->temp);
    rb_mls_init(&seq, RB_PRE_SEQ_POLY);
    for (i = 0; i < RB_PRE_SEQ_LEN; ++i)
        rx->freq[bin(rx, i + RB_PRE_SEQ_OFF)] =
            cf_scale((float)rb_nrz(rb_mls_next(&seq)), rx->freq[bin(rx, i + RB_PRE_SEQ_OFF)]);
    for (i = 0; i < RB_PRE_SEQ_LEN; ++i)
        rx->soft[i] = rb_quantize(RB_BPSK_DIST, 32.0f,
                                  demod_or_erase(rx->freq[bin(rx, i + RB_PRE_SEQ_OFF)],
                                                 rx->freq[bin(rx, i - 1 + RB_PRE_SEQ_OFF)]).re);
    if (!rb_osd_decode(rx->osd, rx->data, rx->soft, rx->generator))
        return ST_FAIL;
    for (i = 0; i < 55; ++i)
        md |= (uint64_t)rb_get_be_bit(rx->data, i) << i;
    for (i = 0; i < 16; ++i)
        cs |= (uint16_t)(rb_get_be_bit(rx->data, i + 55) << i);
    if (rb_crc16_u64(md << 9) != cs)
        return ST_FAIL;
    rx->staged_mode = (int)(md & 255);
    rx->staged_call = md >> 8;
    if (rx->staged_mode && (rx->staged_mode < 14 || rx->staged_mode > 16))
        return ST_NOPE;
    if (rx->staged_call == 0 || rx->staged_call >= 129961739795077ULL) {
        rx->staged_call = 0;
        return ST_NOPE;
    }
    if (!rx->staged_mode)
        return ST_PING;
    return ST_OKAY;
}

/* Decoder::feed for one block (at most extended_length samples, mono). */
static int feed(ribbit_rx *rx, const int16_t *samples, int n)
{
    int i;
    for (i = 0; i < n; ++i) {
        const cf *win = bip_push(rx, analytic(rx, (float)samples[i] / 32768.f));
        if (sc_run(&rx->correlator, win)) {
            rx->stored_cfo_rad = rx->correlator.cfo_rad;
            rx->stored_position = rx->correlator.symbol_pos + rx->accumulated;
            rx->stored_check = 1;
        }
        if (++rx->accumulated == rx->E)
            rx->buf_off = rx->pos0 < rx->pos1 ? rx->pos0 : rx->pos1;
    }
    if (rx->accumulated >= rx->E) {
        rx->accumulated -= rx->E;
        if (rx->stored_check) {
            rx->staged_cfo_rad = rx->stored_cfo_rad;
            rx->staged_position = rx->stored_position;
            rx->staged_check = 1;
            rx->stored_check = 0;
        }
        return 1;
    }
    return 0;
}

static int process(ribbit_rx *rx)
{
    int status = ST_OKAY, i;
    if (rx->staged_check) {
        rx->staged_check = 0;
        status = preamble(rx);
        if (status == ST_OKAY) {
            rx->operation_mode = rx->staged_mode;
            rb_phasor_omega(&rx->osc, -rx->staged_cfo_rad);
            rx->symbol_position = rx->staged_position;
            rx->symbol_number = -1;
            status = ST_SYNC;
        }
    }
    if (rx->symbol_number < RB_SYMBOL_COUNT) {
        const cf *buf = rx->bip + rx->buf_off;
        for (i = 0; i < rx->E; ++i)
            rx->temp[i] = cf_mul(buf[rx->symbol_position + i], rb_phasor_next(&rx->osc));
        rb_fft_run(&rx->fwd, rx->freq, rx->temp);
        if (rx->symbol_number >= 0) {
            for (i = 0; i < RB_PAY_CAR_CNT; ++i)
                rx->cons[i] = demod_or_erase(rx->freq[bin(rx, i + RB_PAY_CAR_OFF)], rx->prev[i]);
            compensate(rx);
            demap(rx);
        }
        if (++rx->symbol_number == RB_SYMBOL_COUNT)
            status = ST_DONE;
        for (i = 0; i < RB_PAY_CAR_CNT; ++i)
            rx->prev[i] = rx->freq[bin(rx, i + RB_PAY_CAR_OFF)];
    }
    return status;
}

static void staged_call(const ribbit_rx *rx, char *call)
{
    static const char alphabet[] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    uint64_t val = rx->staged_call;
    int i;
    for (i = RIBBIT_MAX_CALL - 1; i >= 0; --i, val /= 37)
        call[i] = alphabet[val % 37];
    call[RIBBIT_MAX_CALL] = 0;
}

static ribbit_event *push_event(ribbit_rx *rx, int kind)
{
    ribbit_event *ev;
    if (rx->qcount == QUEUE) {
        rx->qhead = (rx->qhead + 1) % QUEUE;
        --rx->qcount;
        ++rx->dropped;
    }
    ev = &rx->queue[(rx->qhead + rx->qcount) % QUEUE];
    ++rx->qcount;
    memset(ev, 0, sizeof(*ev));
    ev->kind = kind;
    ev->sample = rx->consumed;
    return ev;
}

static void staged_event(ribbit_rx *rx, int kind)
{
    ribbit_event *ev = push_event(rx, kind);
    ev->mode = rx->staged_mode;
    ev->cfo_hz = rx->staged_cfo_rad * ((float)rx->rate / RB_TWO_PI);
    staged_call(rx, ev->call);
}

static void handle(ribbit_rx *rx, int status)
{
    ribbit_event *ev;
    int data_bits, result, i;
    const uint32_t *frozen;
    switch (status) {
    case ST_SYNC:
        staged_event(rx, RIBBIT_EV_SYNC);
        rx->sync_cfo_hz = rx->queue[(rx->qhead + rx->qcount - 1) % QUEUE].cfo_hz;
        staged_call(rx, rx->sync_call);
        break;
    case ST_PING:
        staged_event(rx, RIBBIT_EV_PING);
        break;
    case ST_NOPE:
        staged_event(rx, RIBBIT_EV_NOPE);
        break;
    case ST_DONE:
        frozen = rb_mode_frozen(rx->operation_mode, &data_bits);
        ev = push_event(rx, RIBBIT_EV_FAIL);
        ev->mode = rx->operation_mode;
        ev->cfo_hz = rx->sync_cfo_hz;
        memcpy(ev->call, rx->sync_call, sizeof(ev->call));
        if (!frozen) {
            ev->flips = -1;
            break;
        }
        /* Decoder::fetch */
        result = rb_polar_decode(rx->polar, ev->payload, rx->code, frozen, data_bits);
        {
            rb_xorshift32 scrambler;
            rb_xorshift32_init(&scrambler);
            for (i = 0; i < data_bits / 8; ++i)
                ev->payload[i] ^= (uint8_t)rb_xorshift32_next(&scrambler);
            for (i = data_bits / 8; i < RIBBIT_MAX_PAYLOAD; ++i)
                ev->payload[i] = 0;
        }
        ev->flips = result;
        if (result < 0) {
            memset(ev->payload, 0, sizeof(ev->payload));
            break;
        }
        ev->kind = RIBBIT_EV_DONE;
        while (ev->length < RIBBIT_MAX_PAYLOAD && ev->payload[ev->length])
            ++ev->length;
        break;
    default:
        break;
    }
}

ribbit_rx *ribbit_rx_new(int rate)
{
    ribbit_rx *rx;
    cf *seq;
    rb_mls mls;
    int i, half;
    if (!ribbit_rate_supported(rate))
        return NULL;
    rx = (ribbit_rx *)calloc(1, sizeof(ribbit_rx));
    if (!rx)
        return NULL;
    rx->rate = rate;
    rx->S = (1280 * rate) / 8000;
    rx->G = rx->S / 8;
    rx->E = rx->S + rx->G;
    rx->taps = (((33 * rate) / 8000) & ~3) | 1;
    rx->chunk = rate / 50;
    rx->buffer_length = 4 * rx->E;
    rx->search_position = rx->E;
    rx->hreal = (float *)calloc((size_t)rx->taps, sizeof(float));
    rx->imco = (float *)calloc((size_t)((rx->taps - 1) / 4), sizeof(float));
    rx->bip = (cf *)calloc((size_t)(2 * rx->buffer_length), sizeof(cf));
    rx->temp = (cf *)calloc((size_t)rx->E, sizeof(cf));
    rx->freq = (cf *)calloc((size_t)rx->S, sizeof(cf));
    rx->tse_temp = (float *)calloc((size_t)(((RB_PAY_CAR_CNT - 1) * RB_PAY_CAR_CNT) / 2 + 1),
                                   sizeof(float));
    rx->generator = (int8_t *)calloc(255 * 71, 1);
    rx->pending = (int16_t *)calloc((size_t)rx->chunk, sizeof(int16_t));
    rx->osd = rb_osd_new();
    rx->polar = rb_polar_dec_new();
    if (!rx->hreal || !rx->imco || !rx->bip || !rx->temp || !rx->freq || !rx->tse_temp ||
        !rx->generator || !rx->pending || !rx->osd || !rx->polar ||
        rb_fft_init(&rx->fwd, rx->S, -1)) {
        ribbit_rx_free(rx);
        return NULL;
    }
    /* Decoder(): correlator(corSeq()) */
    half = rx->S / 2;
    seq = rx->freq;
    for (i = 0; i < half; ++i)
        seq[i] = cf_make(0, 0);
    rb_mls_init(&mls, RB_COR_SEQ_POLY);
    for (i = 0; i < RB_COR_SEQ_LEN; ++i)
        seq[(i + RB_COR_SEQ_OFF / 2 + half) % half] = cf_make((float)rb_nrz(rb_mls_next(&mls)), 0);
    if (sc_init(&rx->correlator, rx->search_position, half, rx->G, seq)) {
        ribbit_rx_free(rx);
        return NULL;
    }
    rb_bch_matrix(rx->generator);
    /* Hilbert<cmplx, filter_length>(2) */
    rx->reco = kaiser(2.0f, (rx->taps - 1) / 2, rx->taps);
    for (i = 0; i < (rx->taps - 1) / 4; ++i)
        rx->imco[i] = kaiser(2.0f, (2 * i + 1) + (rx->taps - 1) / 2, rx->taps) * 2 /
                      ((float)(2 * i + 1) * RB_PI);
    /* block_dc.samples(filter_length) */
    rx->dc_a = (float)(rx->taps - 1) / (float)rx->taps;
    rx->dc_b = (1.0f + rx->dc_a) / 2.0f;
    rx->pos0 = 0;
    rx->pos1 = rx->buffer_length;
    rb_phasor_init(&rx->osc);
    rb_phasor_omega_n(&rx->osc, -2000, rate);
    rx->symbol_number = RB_SYMBOL_COUNT;
    rx->symbol_position = rx->search_position + rx->E;
    return rx;
}

void ribbit_rx_free(ribbit_rx *rx)
{
    if (!rx)
        return;
    free(rx->hreal);
    free(rx->imco);
    free(rx->bip);
    free(rx->temp);
    free(rx->freq);
    free(rx->tse_temp);
    free(rx->generator);
    free(rx->pending);
    rb_osd_free(rx->osd);
    rb_polar_dec_free(rx->polar);
    rb_fft_free(&rx->fwd);
    sc_free(&rx->correlator);
    free(rx);
}

int ribbit_rx_process(ribbit_rx *rx, const int16_t *samples, int n)
{
    int i;
    if (!rx || n < 0 || (n && !samples))
        return rx ? rx->qcount : 0;
    for (i = 0; i < n; ++i) {
        rx->pending[rx->npending++] = samples[i];
        ++rx->consumed;
        if (rx->npending == rx->chunk) {
            rx->npending = 0;
            if (feed(rx, rx->pending, rx->chunk))
                handle(rx, process(rx));
        }
    }
    return rx->qcount;
}

int ribbit_rx_next_event(ribbit_rx *rx, ribbit_event *ev)
{
    if (!rx || !rx->qcount)
        return 0;
    if (ev)
        *ev = rx->queue[rx->qhead];
    rx->qhead = (rx->qhead + 1) % QUEUE;
    --rx->qcount;
    return 1;
}

long ribbit_rx_dropped(const ribbit_rx *rx) { return rx ? rx->dropped : 0; }
