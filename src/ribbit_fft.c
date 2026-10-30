/*
 * Port of the aicodix unit circle (unit_circle.hh, Cephes polynomials) and the
 * mixed-radix decimation-in-time FFT (fft.hh). The reference picks radices at compile
 * time through templates; here the same split() rule runs at run time, and every
 * butterfly keeps the reference's operand order, so the results are bit-identical.
 * Radices 2, 3, 4, 5, 7 and 8 cover every size the modem uses (symbol lengths 1280,
 * 2560, 5120, 7056, 7680, their halves and the PAPR sizes).
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <stdlib.h>
#include "ribbit_internal.h"

/* ---- unit_circle.hh ---- */
static const float uc_a1 = 1.0f, uc_a2 = -0.5f;
#define UC(name, val) static const float name = (float)(val)
UC(uc_a3, -1.66666666666666307295E-1);
UC(uc_a4, 4.16666666666665929218E-2);
UC(uc_a5, 8.33333333332211858878E-3);
UC(uc_a6, -1.38888888888730564116E-3);
UC(uc_a7, -1.98412698295895385996E-4);
UC(uc_a8, 2.48015872888517045348E-5);
UC(uc_a9, 2.75573136213857245213E-6);
UC(uc_a10, -2.75573141792967388112E-7);
UC(uc_a11, -2.50507477628578072866E-8);
UC(uc_a12, 2.08757008419747316778E-9);
UC(uc_a13, 1.58962301576546568060E-10);
UC(uc_a14, -1.13585365213876817300E-11);

static float cosine_kernel(float xx)
{
    return uc_a1 + xx * (xx * (xx * (xx * (xx * (xx * (xx * uc_a14 + uc_a12) + uc_a10) +
                                            uc_a8) + uc_a6) + uc_a4) + uc_a2);
}

static float sine_kernel(float x, float xx)
{
    return x * (xx * (xx * (xx * (xx * (xx * (xx * uc_a13 + uc_a11) + uc_a9) + uc_a7) +
                              uc_a5) + uc_a3) + uc_a1);
}

static float cosine_approx(float x) { return cosine_kernel(x * x); }
static float sine_approx(float x) { return sine_kernel(x, x * x); }
static float rad(int n, int N) { return RB_TWO_PI * (float)n / (float)N; }

float rb_uc_cos(int n, int N)
{
    return
        8 * n < -7 * N ? cosine_approx(rad(n + N, N)) :
        8 * n < -5 * N ? -sine_approx(rad(4 * n + 3 * N, 4 * N)) :
        8 * n < -3 * N ? -cosine_approx(rad(2 * n + N, 2 * N)) :
        8 * n < -1 * N ? sine_approx(rad(4 * n + N, 4 * N)) :
        8 * n < 1 * N ? cosine_approx(rad(n, N)) :
        8 * n < 3 * N ? -sine_approx(rad(4 * n - N, 4 * N)) :
        8 * n < 5 * N ? -cosine_approx(rad(2 * n - N, 2 * N)) :
        8 * n < 7 * N ? sine_approx(rad(4 * n - 3 * N, 4 * N)) :
        cosine_approx(rad(n - N, N));
}

float rb_uc_sin(int n, int N)
{
    return
        8 * n < -7 * N ? sine_approx(rad(n + N, N)) :
        8 * n < -5 * N ? cosine_approx(rad(4 * n + 3 * N, 4 * N)) :
        8 * n < -3 * N ? -sine_approx(rad(2 * n + N, 2 * N)) :
        8 * n < -1 * N ? -cosine_approx(rad(4 * n + N, 4 * N)) :
        8 * n < 1 * N ? sine_approx(rad(n, N)) :
        8 * n < 3 * N ? cosine_approx(rad(4 * n - N, 4 * N)) :
        8 * n < 5 * N ? -sine_approx(rad(2 * n - N, 2 * N)) :
        8 * n < 7 * N ? -cosine_approx(rad(4 * n - 3 * N, 4 * N)) :
        sine_approx(rad(n - N, N));
}

/* ---- fft.hh helpers ---- */
static cf cx(int n, int N, cf a) { return cf_scale(rb_uc_cos(n, N), a); }
static cf sx(int n, int N, cf a) { return cf_scale(rb_uc_sin(n, N), a); }

/* rsqrt2<Complex<float>>: Const<Complex<float>>::InvSqrtTwo() is a complex constant, so the
 * product is a full complex multiplication (with a zero imaginary part). */
static cf rsqrt2(cf a) { return cf_mul(cf_make(RB_INV_SQRT_TWO, 0.0f), a); }

static cf fiddle(cf a, cf b)
{
    cf c = cf_add(a, b), d = cf_sub(a, b);
    return cf_make(d.re + c.im, d.im - c.re);
}

static cf twiddle(cf a, cf b) { return cf_make(a.im - b.im, b.re - a.re); }

static int pow2(int N) { return !(N & (N - 1)); }
static int pow4(int N) { return pow2(N) && (N & 0x55555555); }
static int pow8(int N) { return pow2(N) && (N & 0x49249249); }

static int split(int N)
{
    return
        (!(N % 31)) ? 31 :
        (!(N % 29)) ? 29 :
        (!(N % 23)) ? 23 :
        (!(N % 19)) ? 19 :
        (!(N % 17)) ? 17 :
        (!(N % 13)) ? 13 :
        (!(N % 11)) ? 11 :
        (!(N % 7)) ? 7 :
        (!(N % 5)) ? 5 :
        (!(N % 3)) ? 3 :
        (!(N % 8) && pow8(N)) ? 8 :
        (!(N % 8) && pow8(N / 2)) ? 2 :
        (!(N % 4) && pow4(N)) ? 4 :
        (!(N % 2)) ? 2 :
        1;
}

/* The SIGN=+1 butterflies; SIGN=-1 writes the same values to the mirrored outputs
 * (Dit<R, R, STRIDE, TYPE, -1>::dft calls the +1 dft with out1..outR-1 reversed). */
static void dft(int radix, int sign, cf **o, const cf *in)
{
    cf r[8];
    int i;
    switch (radix) {
    case 1:
        r[0] = in[0];
        break;
    case 2:
        r[0] = cf_add(in[0], in[1]);
        r[1] = cf_sub(in[0], in[1]);
        break;
    case 3: {
        cf a1 = cf_add(in[1], in[2]), t1 = twiddle(in[1], in[2]);
        cf c1 = cx(1, 3, a1);
        cf s1 = sx(1, 3, t1);
        r[0] = cf_add(in[0], a1);
        r[1] = cf_sub(cf_add(in[0], c1), s1);
        r[2] = cf_add(cf_add(in[0], c1), s1);
        break;
    }
    case 4: {
        cf a = cf_add(in[0], in[2]), b = cf_sub(in[0], in[2]);
        cf c = cf_add(in[1], in[3]), d = twiddle(in[1], in[3]);
        r[0] = cf_add(a, c);
        r[1] = cf_sub(b, d);
        r[2] = cf_sub(a, c);
        r[3] = cf_add(b, d);
        break;
    }
    case 5: {
        cf a1 = cf_add(in[1], in[4]), t1 = twiddle(in[1], in[4]);
        cf a2 = cf_add(in[2], in[3]), t2 = twiddle(in[2], in[3]);
        cf c1 = cf_add(cx(1, 5, a1), cx(2, 5, a2));
        cf s1 = cf_add(sx(1, 5, t1), sx(2, 5, t2));
        cf c2 = cf_add(cx(2, 5, a1), cx(1, 5, a2));
        cf s2 = cf_sub(sx(2, 5, t1), sx(1, 5, t2));
        r[0] = cf_add(cf_add(in[0], a1), a2);
        r[1] = cf_sub(cf_add(in[0], c1), s1);
        r[2] = cf_sub(cf_add(in[0], c2), s2);
        r[3] = cf_add(cf_add(in[0], c2), s2);
        r[4] = cf_add(cf_add(in[0], c1), s1);
        break;
    }
    case 7: {
        cf a1 = cf_add(in[1], in[6]), t1 = twiddle(in[1], in[6]);
        cf a2 = cf_add(in[2], in[5]), t2 = twiddle(in[2], in[5]);
        cf a3 = cf_add(in[3], in[4]), t3 = twiddle(in[3], in[4]);
        cf c1 = cf_add(cf_add(cx(1, 7, a1), cx(2, 7, a2)), cx(3, 7, a3));
        cf s1 = cf_add(cf_add(sx(1, 7, t1), sx(2, 7, t2)), sx(3, 7, t3));
        cf c2 = cf_add(cf_add(cx(2, 7, a1), cx(3, 7, a2)), cx(1, 7, a3));
        cf s2 = cf_sub(cf_sub(sx(2, 7, t1), sx(3, 7, t2)), sx(1, 7, t3));
        cf c3 = cf_add(cf_add(cx(3, 7, a1), cx(1, 7, a2)), cx(2, 7, a3));
        cf s3 = cf_add(cf_sub(sx(3, 7, t1), sx(1, 7, t2)), sx(2, 7, t3));
        r[0] = cf_add(cf_add(cf_add(in[0], a1), a2), a3);
        r[1] = cf_sub(cf_add(in[0], c1), s1);
        r[2] = cf_sub(cf_add(in[0], c2), s2);
        r[3] = cf_sub(cf_add(in[0], c3), s3);
        r[4] = cf_add(cf_add(in[0], c3), s3);
        r[5] = cf_add(cf_add(in[0], c2), s2);
        r[6] = cf_add(cf_add(in[0], c1), s1);
        break;
    }
    case 8: {
        cf a = cf_add(in[0], in[4]), b = cf_sub(in[0], in[4]);
        cf c = cf_add(in[1], in[5]), d = cf_sub(in[1], in[5]);
        cf e = cf_add(in[2], in[6]), f = twiddle(in[2], in[6]);
        cf g = cf_add(in[3], in[7]), h = cf_sub(in[3], in[7]);
        cf cpg = cf_add(c, g), tcg = twiddle(c, g);
        cf fdh = rsqrt2(fiddle(d, h)), fhd = rsqrt2(fiddle(h, d));
        r[0] = cf_add(cf_add(a, e), cpg);
        r[1] = cf_sub(cf_sub(b, f), fhd);
        r[2] = cf_sub(cf_sub(a, e), tcg);
        r[3] = cf_sub(cf_add(b, f), fdh);
        r[4] = cf_sub(cf_add(a, e), cpg);
        r[5] = cf_add(cf_sub(b, f), fhd);
        r[6] = cf_add(cf_sub(a, e), tcg);
        r[7] = cf_add(cf_add(b, f), fdh);
        break;
    }
    default:
        return;
    }
    *o[0] = r[0];
    for (i = 1; i < radix; ++i)
        *o[sign > 0 ? i : radix - i] = r[i];
}

static void dit(int N, int stride, int sign, cf *out, const cf *in, const cf *z)
{
    int radix = split(N), q, o, i, j, k;
    cf *op[8];
    cf v[8];
    if (radix == N) {
        for (j = 0; j < radix; ++j) {
            op[j] = out + j;
            v[j] = in[j * stride];
        }
        dft(radix, sign, op, v);
        return;
    }
    q = N / radix;
    for (o = 0, i = 0; o < N; o += q, i += stride)
        dit(q, radix * stride, sign, out + o, in + i, z);
    for (k = 0; k < q; ++k) {
        op[0] = out + k;
        v[0] = out[k];
        for (j = 1; j < radix; ++j) {
            op[j] = out + k + j * q;
            v[j] = cf_mul(z[j * stride * k], out[k + j * q]);
        }
        dft(radix, sign, op, v);
    }
}

static int supported(int n)
{
    while (n > 1) {
        int r = split(n);
        if (r != 2 && r != 3 && r != 4 && r != 5 && r != 7 && r != 8)
            return 0;
        n /= r;
    }
    return 1;
}

int rb_fft_init(rb_fft *f, int n, int sign)
{
    int k;
    f->z = NULL;
    if (n < 1 || !supported(n))
        return -1;
    f->n = n;
    f->sign = sign;
    f->z = (cf *)malloc(sizeof(cf) * (size_t)n);
    if (!f->z)
        return -1;
    for (k = 0; k < n; ++k)
        f->z[k] = cf_make(rb_uc_cos(k, n), (float)sign * rb_uc_sin(k, n));
    return 0;
}

void rb_fft_free(rb_fft *f)
{
    free(f->z);
    f->z = NULL;
}

void rb_fft_run(const rb_fft *f, cf *out, const cf *in)
{
    dit(f->n, 1, f->sign, out, in, f->z);
}
