/*
 * Port of the aicodix coding helpers Rattlegram uses: mls.hh (maximum length sequences),
 * crc.hh, bose_chaudhuri_hocquenghem_encoder.hh and osd.hh (BCH(255, 71) preamble code and
 * its ordered statistics decoder, order 2), sort.hh (the stable merge sort OSD uses),
 * polar_encoder.hh and polar_list_decoder.hh with polar_helper.hh (CA-SCL polar code,
 * polar.hh). The list decoder runs 16 paths, the lane count of the reference's SSE4.1 and
 * NEON builds (polar.hh: SIMD<int8_t, 16 / sizeof(int8_t)>), with the NEON lane arithmetic
 * (saturating add and abs, wrapping multiply); the SSE4.1 variant computes the same values
 * for the hard decisions (+1/-1) it is applied to.
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <stdlib.h>
#include <string.h>
#include "ribbit_internal.h"

/* ---- mls.hh ---- */
static int hibit(unsigned n)
{
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return (int)(n ^ (n >> 1));
}

void rb_mls_init(rb_mls *m, int poly)
{
    m->poly = poly;
    m->test = hibit((unsigned)poly) >> 1;
    m->reg = 1;
}

int rb_mls_next(rb_mls *m)
{
    int fb = (m->reg & m->test) != 0;
    m->reg <<= 1;
    m->reg ^= fb * m->poly;
    return fb;
}

/* ---- crc.hh ---- */
uint16_t rb_crc16_u64(uint64_t data)
{
    const uint16_t poly = 0xA8F4;
    uint16_t lut[256], crc = 0;
    int i, j;
    for (j = 0; j < 256; ++j) {
        uint16_t tmp = (uint16_t)j;
        for (i = 8; i; --i)
            tmp = (uint16_t)((tmp >> 1) ^ ((tmp & 1) * poly));
        lut[j] = tmp;
    }
    for (i = 0; i < 8; ++i) {
        uint16_t tmp = (uint16_t)(crc ^ (uint8_t)((data >> (8 * i)) & 255));
        crc = (uint16_t)((crc >> 8) ^ lut[tmp & 255]);
    }
    return crc;
}

void rb_crc32_init(rb_crc32 *c)
{
    const uint32_t poly = 0x8F6E37A0u;
    int i, j;
    for (j = 0; j < 256; ++j) {
        uint32_t tmp = (uint32_t)j;
        for (i = 8; i; --i)
            tmp = (tmp >> 1) ^ ((tmp & 1) * poly);
        c->lut[j] = tmp;
    }
    c->crc = 0;
}

/* ---- BCH(255, 71): minimal polynomials from encoder.hh / decoder.hh ---- */
enum { BCH_N = 255, BCH_K = 71, BCH_NP = BCH_N - BCH_K, BCH_G = ((BCH_NP + 1) + 7) / 8 };
static const int bch_minpolys[24] = {
    0x11D, 0x177, 0x1F3, 0x169, 0x1BD, 0x1E7, 0x12B, 0x1D7,
    0x013, 0x165, 0x18B, 0x163, 0x11B, 0x13F, 0x18D, 0x12D,
    0x15F, 0x1F9, 0x1C3, 0x139, 0x1A9, 0x01F, 0x187, 0x1B1
};

static int degree(int m)
{
    int d = 0;
    while (m >> d)
        ++d;
    return d - 1;
}

/* BoseChaudhuriHocquenghemEncoder<255, 71>::operator() with its constructor. */
void rb_bch_encode(const uint8_t *data, uint8_t *parity)
{
    uint8_t generator[BCH_G];
    int generator_degree = 1, n, i, j, l;
    memset(generator, 0, sizeof(generator));
    rb_set_be_bit(generator, BCH_NP, 1);
    for (n = 0; n < 24; ++n) {
        int m = bch_minpolys[n], m_degree = degree(m);
        for (i = generator_degree; i >= 0; --i) {
            if (!rb_get_be_bit(generator, BCH_NP - i))
                continue;
            rb_set_be_bit(generator, BCH_NP - i, m & 1);
            for (j = 1; j <= m_degree; ++j)
                rb_xor_be_bit(generator, BCH_NP - (i + j), (m >> j) & 1);
        }
        generator_degree += m_degree;
    }
    for (i = 0; i < BCH_NP; ++i)
        rb_set_be_bit(generator, i, rb_get_be_bit(generator, i + 1));
    rb_set_be_bit(generator, BCH_NP, 0);

    for (l = 0; l <= (BCH_NP - 1) / 8; ++l)
        parity[l] = 0;
    for (i = 0; i < BCH_K; ++i) {
        if (rb_get_be_bit(data, i) != rb_get_be_bit(parity, 0)) {
            for (l = 0; l < (BCH_NP - 1) / 8; ++l)
                parity[l] = (uint8_t)(generator[l] ^ ((parity[l] << 1) | (parity[l + 1] >> 7)));
            parity[(BCH_NP - 1) / 8] =
                (uint8_t)(generator[(BCH_NP - 1) / 8] ^ (parity[(BCH_NP - 1) / 8] << 1));
        } else {
            for (l = 0; l < (BCH_NP - 1) / 8; ++l)
                parity[l] = (uint8_t)((parity[l] << 1) | (parity[l + 1] >> 7));
            parity[(BCH_NP - 1) / 8] = (uint8_t)(parity[(BCH_NP - 1) / 8] << 1);
        }
    }
}

/* BoseChaudhuriHocquenghemGenerator<255, 71>::matrix(genmat, true, ...). */
void rb_bch_matrix(int8_t *genmat)
{
    int genpoly_degree = 1, n, i, j, k;
    int8_t *genpoly = genmat;
    for (i = 0; i < BCH_NP; ++i)
        genpoly[i] = 0;
    genpoly[BCH_NP] = 1;
    for (n = 0; n < 24; ++n) {
        int m = bch_minpolys[n], m_degree = degree(m);
        for (i = genpoly_degree; i >= 0; --i) {
            if (!genpoly[BCH_NP - i])
                continue;
            genpoly[BCH_NP - i] = (int8_t)(m & 1);
            for (j = 1; j <= m_degree; ++j)
                genpoly[BCH_NP - (i + j)] ^= (int8_t)((m >> j) & 1);
        }
        genpoly_degree += m_degree;
    }
    for (i = BCH_NP + 1; i < BCH_N; ++i)
        genmat[i] = 0;
    for (j = 1; j < BCH_K; ++j) {
        for (i = 0; i < j; ++i)
            genmat[BCH_N * j + i] = 0;
        for (i = 0; i <= BCH_NP; ++i)
            genmat[(BCH_N + 1) * j + i] = genmat[i];
        for (i = j + BCH_NP + 1; i < BCH_N; ++i)
            genmat[BCH_N * j + i] = 0;
    }
    for (k = BCH_K - 1; k; --k)
        for (j = 0; j < k; ++j)
            if (genmat[BCH_N * j + k])
                for (i = k; i < BCH_N; ++i)
                    genmat[BCH_N * j + i] ^= genmat[BCH_N * k + i];
}

/* ---- osd.hh: OrderedStatisticsDecoder<255, 71, 2> ---- */
enum { OSD_W = (BCH_N + 8 - 1) & ~(8 - 1), MERGE_M = 32 };

struct rb_osd {
    int8_t G[OSD_W * BCH_K];
    int8_t codeword[OSD_W], candidate[OSD_W];
    int8_t softperm[OSD_W];
    int16_t perm[OSD_W];
    int16_t tmp[BCH_N];
};

rb_osd *rb_osd_new(void) { return (rb_osd *)calloc(1, sizeof(rb_osd)); }
void rb_osd_free(rb_osd *o) { free(o); }

/* sort.hh: insertion_sort(a, n, comp) and MergeSort<int16_t, 255, 32> with the OSD's
 * comparator softperm[a] > softperm[b]. */
static int osd_comp(const rb_osd *o, int a, int b) { return o->softperm[a] > o->softperm[b]; }

static void osd_insertion(const rb_osd *o, int16_t *a, int n)
{
    int i, j;
    for (i = 1; i < n; ++i) {
        int16_t t = a[i];
        for (j = i; j > 0 && osd_comp(o, t, a[j - 1]); --j)
            a[j] = a[j - 1];
        a[j] = t;
    }
}

static void osd_sort(rb_osd *o, int16_t *a, int n)
{
    int i, l;
    for (i = 0; i < n; i += MERGE_M)
        osd_insertion(o, a + i, i > n - MERGE_M ? n - i : MERGE_M);
    for (l = MERGE_M; l < n; l *= 2) {
        for (i = 0; i < n; i += 2 * l) {
            int left = i, right = i + l, end = i + 2 * l, x, y, k;
            if (right > n)
                right = n;
            if (end > n)
                end = n;
            for (x = left, y = right, k = left; k < end; ++k)
                o->tmp[k] = (x >= right || (y < end && osd_comp(o, a[y], a[x]))) ? a[y++] : a[x++];
        }
        for (i = 0; i < n; ++i)
            a[i] = o->tmp[i];
    }
}

static void swap8(int8_t *a, int8_t *b) { int8_t t = *a; *a = *b; *b = t; }

static void osd_row_echelon(rb_osd *o)
{
    int8_t *G = o->G;
    int k, j, i, h;
    for (k = 0; k < BCH_K; ++k) {
        for (j = k; j < BCH_K; ++j) {
            if (G[OSD_W * j + k]) {
                for (i = k; j != k && i < BCH_N; ++i)
                    swap8(&G[OSD_W * j + i], &G[OSD_W * k + i]);
                break;
            }
        }
        for (j = k + 1; !G[OSD_W * k + k] && j < BCH_N; ++j) {
            for (h = k; h < BCH_K; ++h) {
                if (G[OSD_W * h + j]) {
                    int16_t t = o->perm[k];
                    o->perm[k] = o->perm[j];
                    o->perm[j] = t;
                    for (i = 0; i < BCH_K; ++i)
                        swap8(&G[OSD_W * i + k], &G[OSD_W * i + j]);
                    for (i = k; h != k && i < BCH_N; ++i)
                        swap8(&G[OSD_W * h + i], &G[OSD_W * k + i]);
                    break;
                }
            }
        }
        for (j = k + 1; j < BCH_K; ++j)
            if (G[OSD_W * j + k])
                for (i = k; i < BCH_N; ++i)
                    G[OSD_W * j + i] ^= G[OSD_W * k + i];
    }
}

static void osd_systematic(rb_osd *o)
{
    int8_t *G = o->G;
    int k, j, i;
    for (k = BCH_K - 1; k; --k)
        for (j = 0; j < k; ++j)
            if (G[OSD_W * j + k])
                for (i = k; i < BCH_N; ++i)
                    G[OSD_W * j + i] ^= G[OSD_W * k + i];
}

static void osd_encode(rb_osd *o)
{
    int i, j;
    for (i = BCH_K; i < BCH_N; ++i)
        o->codeword[i] = (int8_t)(o->codeword[0] & o->G[i]);
    for (j = 1; j < BCH_K; ++j)
        for (i = BCH_K; i < BCH_N; ++i)
            o->codeword[i] ^= (int8_t)(o->codeword[j] & o->G[OSD_W * j + i]);
}

static void osd_flip(rb_osd *o, int j)
{
    int i;
    for (i = 0; i < OSD_W; ++i)
        o->codeword[i] ^= o->G[OSD_W * j + i];
}

static int osd_metric(const rb_osd *o)
{
    int sum = 0, i;
    for (i = 0; i < OSD_W; ++i)
        sum += (1 - 2 * o->codeword[i]) * o->softperm[i];
    return sum;
}

static void osd_update(rb_osd *o, int *best, int *next)
{
    int met = osd_metric(o);
    if (met > *best) {
        *next = *best;
        *best = met;
        memcpy(o->candidate, o->codeword, BCH_N);
    } else if (met > *next) {
        *next = met;
    }
}

static int8_t max127(int8_t v) { return v < -127 ? (int8_t)-127 : v; }

int rb_osd_decode(rb_osd *o, uint8_t *hard, const int8_t *soft, const int8_t *genmat)
{
    int i, j, a, b, best, next = -1;
    for (i = 0; i < BCH_N; ++i)
        o->perm[i] = (int16_t)i;
    for (i = 0; i < BCH_N; ++i) {
        int8_t v = max127(soft[i]);
        o->softperm[i] = (int8_t)(v < 0 ? -v : v);
    }
    osd_sort(o, o->perm, BCH_N);
    for (j = 0; j < BCH_K; ++j)
        for (i = 0; i < BCH_N; ++i)
            o->G[OSD_W * j + i] = genmat[BCH_N * j + o->perm[i]];
    osd_row_echelon(o);
    osd_systematic(o);
    for (i = 0; i < BCH_N; ++i)
        o->softperm[i] = max127(soft[o->perm[i]]);
    for (i = BCH_N; i < OSD_W; ++i)
        o->softperm[i] = 0;
    for (i = 0; i < BCH_K; ++i)
        o->codeword[i] = (int8_t)(o->softperm[i] < 0);
    osd_encode(o);
    memcpy(o->candidate, o->codeword, BCH_N);
    best = osd_metric(o);
    for (a = 0; a < BCH_K; ++a) {
        osd_flip(o, a);
        osd_update(o, &best, &next);
        for (b = a + 1; b < BCH_K; ++b) {
            osd_flip(o, b);
            osd_update(o, &best, &next);
            osd_flip(o, b);
        }
        osd_flip(o, a);
    }
    for (i = 0; i < BCH_N; ++i)
        rb_set_be_bit(hard, o->perm[i], o->candidate[i]);
    return best != next;
}

/* ---- polar_encoder.hh: PolarSysEnc<int8_t> / polar.hh: PolarEncoder ---- */
static int frozen_bit(const uint32_t *bits, int idx) { return (bits[idx / 32] >> (idx % 32)) & 1; }

void rb_polar_encode(int8_t *code, const uint8_t *message, const uint32_t *frozen, int data_bits)
{
    int8_t mesg[1360 + 32];
    const int8_t *msg = mesg;
    rb_crc32 crc;
    int i, h, j, round;
    for (i = 0; i < data_bits; ++i)
        mesg[i] = (int8_t)rb_nrz(rb_get_le_bit(message, i));
    rb_crc32_init(&crc);
    for (i = 0; i < data_bits / 8; ++i)
        rb_crc32_byte(&crc, message[i]);
    for (i = 0; i < 32; ++i)
        mesg[i + data_bits] = (int8_t)rb_nrz((crc.crc >> i) & 1);
    for (round = 0; round < 2; ++round) {
        for (i = 0; i < RB_CODE_LEN; i += 2) {
            int8_t msg0, msg1;
            if (round == 0) {
                msg0 = frozen_bit(frozen, i) ? 1 : *msg++;
                msg1 = frozen_bit(frozen, i + 1) ? 1 : *msg++;
            } else {
                msg0 = frozen_bit(frozen, i) ? 1 : code[i];
                msg1 = frozen_bit(frozen, i + 1) ? 1 : code[i + 1];
            }
            code[i] = (int8_t)(msg0 * msg1);
            code[i + 1] = msg1;
        }
        for (h = 2; h < RB_CODE_LEN; h *= 2)
            for (i = 0; i < RB_CODE_LEN; i += 2 * h)
                for (j = i; j < i + h; ++j)
                    code[j] = (int8_t)(code[j] * code[j + h]);
    }
}

/* ---- polar_list_decoder.hh: PolarListDecoder<SIMD<int8_t, 16>, 11> ---- */
enum { L = 16 };
typedef struct { int8_t v[L]; } lane;
typedef struct { uint8_t v[L]; } lmap;

struct rb_polar_dec {
    lane soft[2 * RB_CODE_LEN];
    lane hard[RB_CODE_LEN];
    lmap maps[RB_CODE_LEN];
    lane mesg[1360 + 32];
    lane mess[RB_CODE_LEN];
    int metric[L];
    int count;
};

rb_polar_dec *rb_polar_dec_new(void) { return (rb_polar_dec *)calloc(1, sizeof(rb_polar_dec)); }
void rb_polar_dec_free(rb_polar_dec *d) { free(d); }

static int8_t sat8(int v) { return (int8_t)(v < -128 ? -128 : v > 127 ? 127 : v); }
static int8_t qabs8(int8_t a) { return a == -128 ? (int8_t)127 : (int8_t)(a < 0 ? -a : a); }
static int8_t signum8(int8_t a) { return (int8_t)((a > 0) - (a < 0)); }

static lane l_dup(int8_t a)
{
    lane t;
    int k;
    for (k = 0; k < L; ++k)
        t.v[k] = a;
    return t;
}

static lane l_one(void) { return l_dup(1); }

static lane l_prod(lane a, lane b)
{
    lane t;
    int k;
    for (k = 0; k < L; ++k) {
        int8_t ma = qabs8(a.v[k]), mb = qabs8(b.v[k]);
        t.v[k] = (int8_t)(signum8(a.v[k]) * signum8(b.v[k]) * (ma < mb ? ma : mb));
    }
    return t;
}

static lane l_madd(lane a, lane b, lane c)
{
    lane t;
    int k;
    for (k = 0; k < L; ++k) {
        int8_t bb = b.v[k] < -127 ? (int8_t)-127 : b.v[k];
        int8_t s = sat8((int8_t)(a.v[k] * bb) + c.v[k]);
        t.v[k] = s < -127 ? (int8_t)-127 : s;
    }
    return t;
}

static lane l_qmul(lane a, lane b)
{
    lane t;
    int k;
    for (k = 0; k < L; ++k)
        t.v[k] = (int8_t)(a.v[k] * b.v[k]);
    return t;
}

static lane l_shuf(lane a, lmap m)
{
    lane t;
    int k;
    for (k = 0; k < L; ++k)
        t.v[k] = m.v[k] < L ? a.v[m.v[k]] : 0;
    return t;
}

static lmap m_shuf(lmap a, lmap m)
{
    lmap t;
    int k;
    for (k = 0; k < L; ++k)
        t.v[k] = m.v[k] < L ? a.v[m.v[k]] : 0;
    return t;
}

static lmap m_ident(void)
{
    lmap t;
    int k;
    for (k = 0; k < L; ++k)
        t.v[k] = (uint8_t)k;
    return t;
}

/* PolarListNode<TYPE, M>::rate0 for a node of n = 1 << M bits. */
static lmap node_rate0(rb_polar_dec *d, int n, lane *hard, const lane *soft)
{
    int i, k;
    for (i = 0; i < n; ++i)
        hard[i] = l_one();
    for (i = 0; i < n; ++i)
        for (k = 0; k < L; ++k)
            if (soft[i + n].v[k] < 0)
                d->metric[k] -= soft[i + n].v[k];
    return m_ident();
}

/* PolarListNode<TYPE, 0>::rate1. */
static lmap node_rate1(rb_polar_dec *d, lane *message, lane *hard, const lane *soft)
{
    lane sft = soft[1], hrd;
    int fork[2 * L], perm[2 * L];
    lmap map;
    int i, j, k;
    for (k = 0; k < L; ++k)
        fork[2 * k] = fork[2 * k + 1] = d->metric[k];
    for (k = 0; k < L; ++k)
        if (sft.v[k] < 0)
            fork[2 * k] -= sft.v[k];
        else
            fork[2 * k + 1] += sft.v[k];
    /* sort.hh insertion_sort(INDEX *p, TYPE *a, n) */
    perm[0] = 0;
    for (i = 1; i < 2 * L; ++i) {
        int t = fork[i];
        for (j = i; j > 0 && t < fork[j - 1]; --j) {
            fork[j] = fork[j - 1];
            perm[j] = perm[j - 1];
        }
        fork[j] = t;
        perm[j] = i;
    }
    for (k = 0; k < L; ++k)
        d->metric[k] = fork[k];
    for (k = 0; k < L; ++k)
        map.v[k] = (uint8_t)(perm[k] >> 1);
    for (k = 0; k < L; ++k)
        hrd.v[k] = (int8_t)(1 - 2 * (perm[k] & 1));
    message[d->count] = hrd;
    d->maps[d->count] = map;
    ++d->count;
    *hard = hrd;
    return map;
}

/* PolarListTree<TYPE, M>::decode. For M > 6 the frozen set is the word array fz; M == 6
 * reads fz[0], fz[1]; M <= 5 the 32-bit word fv (polar_list_decoder.hh). */
static lmap tree(rb_polar_dec *d, int M, lane *message, lane *hard, lane *soft,
                 const uint32_t *fz, uint32_t fv)
{
    int N = 1 << M, i;
    lmap lm, rm;
    if (M == 1) {
        soft[1] = l_prod(soft[2], soft[3]);
        lm = (fv & 1) ? node_rate0(d, 1, hard, soft) : node_rate1(d, message, hard, soft);
        soft[1] = l_madd(hard[0], l_shuf(soft[2], lm), l_shuf(soft[3], lm));
        rm = (fv >> 1) ? node_rate0(d, 1, hard + 1, soft)
                       : node_rate1(d, message, hard + 1, soft);
        hard[0] = l_qmul(l_shuf(hard[0], rm), hard[1]);
        return m_shuf(lm, rm);
    }
    for (i = 0; i < N / 2; ++i)
        soft[i + N / 2] = l_prod(soft[i + N], soft[i + N / 2 + N]);
    if (M > 6) {
        lm = tree(d, M - 1, message, hard, soft, fz, 0);
    } else if (M == 6) {
        lm = fz[0] == 0xffffffffu ? node_rate0(d, N / 2, hard, soft)
                                  : tree(d, M - 1, message, hard, soft, NULL, fz[0]);
    } else {
        uint32_t mask = (1u << (1 << (M - 1))) - 1;
        lm = (fv & mask) == mask ? node_rate0(d, N / 2, hard, soft)
                                 : tree(d, M - 1, message, hard, soft, NULL, fv & mask);
    }
    for (i = 0; i < N / 2; ++i)
        soft[i + N / 2] = l_madd(hard[i], l_shuf(soft[i + N], lm), l_shuf(soft[i + N / 2 + N], lm));
    if (M > 6) {
        rm = tree(d, M - 1, message, hard + N / 2, soft, fz + N / 2 / 32, 0);
    } else if (M == 6) {
        rm = fz[1] == 0xffffffffu ? node_rate0(d, N / 2, hard + N / 2, soft)
                                  : tree(d, M - 1, message, hard + N / 2, soft, NULL, fz[1]);
    } else {
        uint32_t mask = (1u << (1 << (M - 1))) - 1;
        rm = (fv >> (N / 2)) == mask ? node_rate0(d, N / 2, hard + N / 2, soft)
                                     : tree(d, M - 1, message, hard + N / 2, soft, NULL,
                                            fv >> (N / 2));
    }
    for (i = 0; i < N / 2; ++i)
        hard[i] = l_qmul(l_shuf(hard[i], rm), hard[i + N / 2]);
    return m_shuf(lm, rm);
}

static void list_decode(rb_polar_dec *d, const int8_t *codeword, const uint32_t *frozen)
{
    int i, k;
    int length = RB_CODE_LEN;
    lmap acc;
    d->count = 0;
    d->metric[0] = 0;
    for (k = 1; k < L; ++k)
        d->metric[k] = 1000000;
    for (i = 0; i < length; ++i)
        d->soft[length + i] = l_dup(codeword[i]);
    tree(d, RB_CODE_ORDER, d->mesg, d->hard, d->soft, frozen, 0);
    acc = d->maps[d->count - 1];
    for (i = d->count - 2; i >= 0; --i) {
        d->mesg[i] = l_shuf(d->mesg[i], acc);
        acc = m_shuf(d->maps[i], acc);
    }
}

int rb_polar_decode(rb_polar_dec *d, uint8_t *message, const int8_t *code,
                    const uint32_t *frozen, int data_bits)
{
    int crc_bits = data_bits + 32, best = -1, flips = 0;
    int i, j, k, h;
    const lane *msg;
    rb_crc32 crc;
    list_decode(d, code, frozen);
    /* systematic(): CODE::PolarEncoder<mesg_type> (non-systematic) of the decoded paths */
    msg = d->mesg;
    for (i = 0; i < RB_CODE_LEN; i += 2) {
        lane msg0 = frozen_bit(frozen, i) ? l_one() : *msg++;
        lane msg1 = frozen_bit(frozen, i + 1) ? l_one() : *msg++;
        d->mess[i] = l_qmul(msg0, msg1);
        d->mess[i + 1] = msg1;
    }
    for (h = 2; h < RB_CODE_LEN; h *= 2)
        for (i = 0; i < RB_CODE_LEN; i += 2 * h)
            for (j = i; j < i + h; ++j)
                d->mess[j] = l_qmul(d->mess[j], d->mess[j + h]);
    for (i = 0, j = 0; i < RB_CODE_LEN && j < crc_bits; ++i)
        if (!frozen_bit(frozen, i))
            d->mesg[j++] = d->mess[i];
    rb_crc32_init(&crc);
    for (k = 0; k < L; ++k) {
        crc.crc = 0;
        for (i = 0; i < crc_bits; ++i)
            rb_crc32_bit(&crc, d->mesg[i].v[k] < 0);
        if (crc.crc == 0) {
            best = k;
            break;
        }
    }
    if (best < 0)
        return -1;
    for (i = 0, j = 0; i < data_bits; ++i, ++j) {
        int received, decoded;
        while (frozen_bit(frozen, j))
            ++j;
        received = code[j] < 0;
        decoded = d->mesg[i].v[best] < 0;
        flips += received != decoded;
        rb_set_le_bit(message, i, decoded);
    }
    return flips;
}
