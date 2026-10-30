/*
 * ribbit_tool -- command line front end of ribbit-core with the same interface and output
 * format as tools/ref_tool.cpp (the C++ reference), so the two can be diffed.
 *
 *   ribbit_tool tx RATE CARRIER NOISE FANCY CALL PAYLOAD_HEX OUT.raw   (PAYLOAD_HEX "-": ping)
 *   ribbit_tool rx RATE IN.raw                                          (events on stdout)
 *
 * Raw files are mono signed 16-bit little-endian samples.
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ribbit.h"

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
    return hex[0] ? -1 : n;
}

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

static int do_tx(int argc, char **argv)
{
    unsigned char payload[RIBBIT_MAX_PAYLOAD];
    int rate, len, n;
    ribbit_tx *tx;
    int16_t *buf;
    FILE *f;
    if (argc != 9)
        return 2;
    rate = atoi(argv[2]);
    len = unhex(argv[7], payload, RIBBIT_MAX_PAYLOAD);
    tx = ribbit_tx_new(rate);
    if (!tx || len < 0 ||
        ribbit_tx_configure(tx, payload, len, argv[6], atoi(argv[3]), atoi(argv[4]),
                            atoi(argv[5]))) {
        fprintf(stderr, "bad arguments\n");
        return 1;
    }
    n = ribbit_symbol_samples(rate);
    buf = (int16_t *)malloc(sizeof(int16_t) * (size_t)n);
    f = fopen(argv[8], "wb");
    if (!buf || !f)
        return 1;
    while (ribbit_tx_produce(tx, buf))
        fwrite(buf, sizeof(int16_t), (size_t)n, f);
    fclose(f);
    free(buf);
    ribbit_tx_free(tx);
    return 0;
}

static int do_rx(int argc, char **argv)
{
    ribbit_rx *rx;
    ribbit_event ev;
    int16_t buf[4096];
    size_t n;
    FILE *f;
    if (argc != 4)
        return 2;
    rx = ribbit_rx_new(atoi(argv[2]));
    f = fopen(argv[3], "rb");
    if (!rx || !f) {
        fprintf(stderr, "bad arguments\n");
        return 1;
    }
    while ((n = fread(buf, sizeof(int16_t), sizeof(buf) / sizeof(buf[0]), f)) > 0) {
        ribbit_rx_process(rx, buf, (int)n);
        while (ribbit_rx_next_event(rx, &ev)) {
            int i;
            printf("%s mode=%d call=[%s] cfo=%a flips=%d len=%d sample=%lld payload=",
                   kind_name(ev.kind), ev.mode, ev.call, (double)ev.cfo_hz, ev.flips,
                   ev.length, ev.sample);
            for (i = 0; i < ev.length; ++i)
                printf("%02x", ev.payload[i]);
            printf("\n");
        }
    }
    fclose(f);
    ribbit_rx_free(rx);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && !strcmp(argv[1], "tx"))
        return do_tx(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "rx"))
        return do_rx(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "version")) {
        printf("%s\n", ribbit_version());
        return 0;
    }
    fprintf(stderr, "usage: ribbit_tool tx RATE CARRIER NOISE FANCY CALL PAYLOAD_HEX OUT.raw\n"
                    "       ribbit_tool rx RATE IN.raw\n");
    return 2;
}
