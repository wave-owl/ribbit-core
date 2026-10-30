/*
 * ribbit-core tests: goldens from the C++ reference (tests/golden, written by
 * tools/ref_check.py --goldens), round trips over every rate and mode, chunking
 * independence and argument checks.
 *
 *   test_ribbit tests/golden
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ribbit.h"

static int failures, checks;

#define CHECK(cond, ...) do { ++checks; if (!(cond)) { ++failures; \
    fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); } } while (0)

static const char *kind_name(int kind)
{
    switch (kind) {
    case RIBBIT_EV_SYNC: return "SYNC";
    case RIBBIT_EV_DONE: return "DONE";
    case RIBBIT_EV_FAIL: return "FAIL";
    case RIBBIT_EV_PING: return "PING";
    case RIBBIT_EV_NOPE: return "NOPE";
    default: return "?";
    }
}

/* Appends the event in tools/ref_tool.cpp's line format. */
static void format_event(char *out, size_t cap, const ribbit_event *ev)
{
    size_t n = strlen(out);
    int i;
    n += (size_t)snprintf(out + n, cap - n, "%s mode=%d call=[%s] cfo=%a flips=%d len=%d "
                          "sample=%lld payload=", kind_name(ev->kind), ev->mode, ev->call,
                          (double)ev->cfo_hz, ev->flips, ev->length, ev->sample);
    for (i = 0; i < ev->length && n + 3 < cap; ++i)
        n += (size_t)snprintf(out + n, cap - n, "%02x", ev->payload[i]);
    snprintf(out + n, cap - n, "\n");
}

static int16_t *read_raw(const char *path, long *count)
{
    FILE *f = fopen(path, "rb");
    long size;
    int16_t *buf;
    if (!f)
        return NULL;
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    buf = (int16_t *)malloc((size_t)size + 2);
    if (buf && fread(buf, 1, (size_t)size, f) != (size_t)size) {
        free(buf);
        buf = NULL;
    }
    fclose(f);
    *count = size / 2;
    return buf;
}

static char *read_text(const char *path)
{
    FILE *f = fopen(path, "rb");
    long size;
    char *buf;
    if (!f)
        return NULL;
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    buf = (char *)calloc((size_t)size + 1, 1);
    if (buf && fread(buf, 1, (size_t)size, f) != (size_t)size) {
        free(buf);
        buf = NULL;
    }
    fclose(f);
    return buf;
}

/* Decodes samples in blocks of `block` (0: all at once) into the event text. */
static char *decode(int rate, const int16_t *s, long n, long block)
{
    size_t cap = 65536;
    char *out = (char *)calloc(cap, 1);
    ribbit_rx *rx = ribbit_rx_new(rate);
    ribbit_event ev;
    long pos = 0;
    if (!out || !rx)
        return out;
    if (block <= 0)
        block = n;
    while (pos < n) {
        long m = n - pos < block ? n - pos : block;
        ribbit_rx_process(rx, s + pos, (int)m);
        pos += m;
        while (ribbit_rx_next_event(rx, &ev))
            format_event(out, cap, &ev);
    }
    ribbit_rx_free(rx);
    return out;
}

static int16_t *render(int rate, const unsigned char *payload, int len, const char *call,
                       int carrier, int noise, int fancy, long *count)
{
    ribbit_tx *tx = ribbit_tx_new(rate);
    int sym = ribbit_symbol_samples(rate);
    long cap = 32L * sym, n = 0;
    int16_t *buf = (int16_t *)calloc((size_t)cap, sizeof(int16_t));
    *count = 0;
    if (!tx || !buf || ribbit_tx_configure(tx, payload, len, call, carrier, noise, fancy)) {
        ribbit_tx_free(tx);
        free(buf);
        return NULL;
    }
    while (n + sym <= cap && ribbit_tx_produce(tx, buf + n))
        n += sym;
    ribbit_tx_free(tx);
    *count = n;
    return buf;
}

static int unhex(const char *hex, unsigned char *out, int max)
{
    int n = 0;
    if (!strcmp(hex, "-"))
        return 0;
    while (hex[0] && hex[1] && n < max) {
        unsigned v;
        if (sscanf(hex, "%2x", &v) != 1)
            return -1;
        out[n++] = (unsigned char)v;
        hex += 2;
    }
    return n;
}

static void test_goldens(const char *dir)
{
    char path[1024], line[2048];
    FILE *idx;
    int rows = 0;
    snprintf(path, sizeof(path), "%s/index.tsv", dir);
    idx = fopen(path, "r");
    CHECK(idx != NULL, "missing %s", path);
    if (!idx)
        return;
    while (fgets(line, sizeof(line), idx)) {
        char name[256], kind[8], call[32] = "", hex[1024] = "-";
        int rate, carrier = 0, noise = 0, fancy = 0;
        long n_raw = 0, n_ch = 0;
        int16_t *raw, *ch;
        char *want, *got;
        if (line[0] == '#')
            continue;
        if (sscanf(line, "%255s %d %7s", name, &rate, kind) != 3)
            continue;
        ++rows;
        snprintf(path, sizeof(path), "%s/%s.raw", dir, name);
        raw = read_raw(path, &n_raw);
        CHECK(raw != NULL, "missing %s", path);
        if (!raw)
            continue;
        if (!strcmp(kind, "tx")) {
            unsigned char payload[RIBBIT_MAX_PAYLOAD];
            long n_tx;
            int16_t *mine;
            int len;
            char *p = strchr(line, '\t');
            p = p ? strchr(p + 1, '\t') : NULL;
            p = p ? strchr(p + 1, '\t') : NULL;
            CHECK(p && sscanf(p + 1, "%31s %d %d %d %1023s", call, &carrier, &noise, &fancy,
                              hex) == 5, "bad index row %s", name);
            len = unhex(hex, payload, RIBBIT_MAX_PAYLOAD);
            mine = render(rate, payload, len, call, carrier, noise, fancy, &n_tx);
            CHECK(mine && n_tx == n_raw && !memcmp(mine, raw, sizeof(int16_t) * (size_t)n_raw),
                  "TX %s differs from the reference (%ld vs %ld samples)", name, n_tx, n_raw);
            free(mine);
            snprintf(path, sizeof(path), "%s/%s.ch.raw", dir, name);
            ch = read_raw(path, &n_ch);
            CHECK(ch != NULL, "missing %s", path);
        } else {
            ch = raw;
            raw = NULL;
            n_ch = n_raw;
        }
        snprintf(path, sizeof(path), "%s/%s.events", dir, name);
        want = read_text(path);
        got = ch ? decode(rate, ch, n_ch, 0) : NULL;
        CHECK(want && got && !strcmp(want, got), "RX %s differs from the reference:\n"
              "  want: %s  got:  %s", name, want ? want : "?", got ? got : "?");
        free(want);
        free(got);
        free(raw);
        free(ch);
    }
    fclose(idx);
    CHECK(rows >= 10, "only %d golden rows", rows);
}

static void test_round_trips(void)
{
    static const int rates[] = {8000, 16000, 32000, 44100, 48000};
    static const int lengths[] = {0, 1, 85, 86, 128, 129, 170};
    unsigned char payload[RIBBIT_MAX_PAYLOAD];
    size_t r, l;
    int i;
    for (i = 0; i < RIBBIT_MAX_PAYLOAD; ++i)
        payload[i] = (unsigned char)(0x20 + (i * 13) % 95);
    for (r = 0; r < sizeof(rates) / sizeof(rates[0]); ++r) {
        for (l = 0; l < sizeof(lengths) / sizeof(lengths[0]); ++l) {
            int rate = rates[r], len = lengths[l], sym = ribbit_symbol_samples(rate);
            long n = 0, total;
            int16_t *sig = render(rate, payload, len, "DG6FL", 1500, 0, 0, &n);
            int16_t *ch;
            ribbit_rx *rx = ribbit_rx_new(rate);
            ribbit_event ev;
            int done = 0, synced = 0, pinged = 0;
            int want_mode = len == 0 ? 0 : len <= 85 ? 16 : len <= 128 ? 15 : 14;
            CHECK(sig && n > 0, "render %d/%d", rate, len);
            total = n + 6L * sym;
            ch = (int16_t *)calloc((size_t)total, sizeof(int16_t));
            if (!sig || !ch || !rx) {
                free(sig);
                free(ch);
                ribbit_rx_free(rx);
                continue;
            }
            memcpy(ch + 3L * sym, sig, sizeof(int16_t) * (size_t)n);
            ribbit_rx_process(rx, ch, (int)total);
            while (ribbit_rx_next_event(rx, &ev)) {
                if (ev.kind == RIBBIT_EV_SYNC)
                    synced = ev.mode == want_mode && !strcmp(ev.call, "    DG6FL");
                if (ev.kind == RIBBIT_EV_PING)
                    pinged = !strcmp(ev.call, "    DG6FL");
                if (ev.kind == RIBBIT_EV_DONE)
                    done = ev.length == len && !memcmp(ev.payload, payload, (size_t)len) &&
                           ev.cfo_hz > 1499.0f && ev.cfo_hz < 1501.0f;
            }
            if (len)
                CHECK(synced && done, "round trip rate %d len %d", rate, len);
            else
                CHECK(pinged, "ping rate %d", rate);
            ribbit_rx_free(rx);
            free(sig);
            free(ch);
        }
    }
}

static void test_chunking(void)
{
    long n = 0;
    int16_t *sig = render(8000, (const unsigned char *)"chunk test", 10, "DG6FL", 1500, 1, 1, &n);
    int16_t *ch = (int16_t *)calloc((size_t)(n + 16000), sizeof(int16_t));
    char *a, *b, *c;
    if (!sig || !ch) {
        CHECK(0, "render");
        return;
    }
    memcpy(ch + 4000, sig, sizeof(int16_t) * (size_t)n);
    a = decode(8000, ch, n + 16000, 0);
    b = decode(8000, ch, n + 16000, 1);
    c = decode(8000, ch, n + 16000, 1037);
    CHECK(a && b && c && !strcmp(a, b) && !strcmp(a, c) && strstr(a, "DONE"),
          "events depend on the caller's block size:\n%s---\n%s---\n%s", a, b, c);
    free(a);
    free(b);
    free(c);
    free(sig);
    free(ch);
}

static void test_arguments(void)
{
    ribbit_tx *tx = ribbit_tx_new(8000);
    unsigned char nul[3] = {'a', 0, 'b'};
    unsigned char big[RIBBIT_MAX_PAYLOAD + 1];
    int16_t out[1440];
    memset(big, 'x', sizeof(big));
    CHECK(ribbit_rx_new(11025) == NULL, "11025 Hz accepted");
    CHECK(ribbit_tx_new(12000) == NULL, "12000 Hz accepted");
    CHECK(ribbit_symbol_samples(8000) == 1440, "symbol samples 8000");
    CHECK(ribbit_symbol_samples(48000) == 8640, "symbol samples 48000");
    CHECK(ribbit_symbol_samples(12000) == 0, "symbol samples 12000");
    CHECK(tx != NULL, "tx 8000");
    if (!tx)
        return;
    CHECK(ribbit_tx_configure(tx, big, RIBBIT_MAX_PAYLOAD + 1, "DG6FL", 1500, 0, 0) == -1,
          "171 bytes accepted");
    CHECK(ribbit_tx_configure(tx, nul, 3, "DG6FL", 1500, 0, 0) == -1, "NUL accepted");
    CHECK(ribbit_tx_configure(tx, big, 5, "DG6FL/P", 1500, 0, 0) == -1, "'/' accepted");
    CHECK(ribbit_tx_configure(tx, big, 5, "ABCDEFGHIJ", 1500, 0, 0) == -1, "10 chars accepted");
    CHECK(ribbit_tx_configure(tx, big, 5, "", 1500, 0, 0) == -1, "empty call accepted");
    CHECK(ribbit_tx_configure(tx, big, 5, "dg6fl", 1500, 0, 0) == 0, "lower case refused");
    CHECK(ribbit_tx_produce(tx, out) == 1, "produce");
    CHECK(ribbit_rx_next_event(NULL, NULL) == 0, "NULL rx");
    CHECK(ribbit_sizeof(RIBBIT_SIZEOF_EVENT) == sizeof(ribbit_event), "sizeof event");
    CHECK(ribbit_sizeof(99) == 0, "sizeof unknown");
    CHECK(ribbit_abi() == RIBBIT_ABI, "abi");
    CHECK(strstr(ribbit_version(), "contract=ribbit-c1 abi=1 tx=1") != NULL, "version %s",
          ribbit_version());
    ribbit_tx_free(tx);
}

int main(int argc, char **argv)
{
    test_arguments();
    test_round_trips();
    test_chunking();
    if (argc > 1)
        test_goldens(argv[1]);
    printf("%s: %d checks, %d failures\n", ribbit_version(), checks, failures);
    return failures ? 1 : 0;
}
