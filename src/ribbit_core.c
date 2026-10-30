/*
 * ribbit-core: version, contract and the tables shared by receiver and transmitter.
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include "ribbit.h"
#include "ribbit_internal.h"

const char *ribbit_version(void)
{
    return "ribbit-core " RIBBIT_VERSION_STRING " contract=ribbit-c1 abi=1 tx=1 list=16";
}

int ribbit_abi(void) { return RIBBIT_ABI; }

size_t ribbit_sizeof(int what)
{
    switch (what) {
    case RIBBIT_SIZEOF_EVENT:
        return sizeof(ribbit_event);
    default:
        return 0;
    }
}

int ribbit_rate_supported(int rate)
{
    return rate == 8000 || rate == 16000 || rate == 32000 || rate == 44100 || rate == 48000;
}

int ribbit_symbol_samples(int rate)
{
    int s;
    if (!ribbit_rate_supported(rate))
        return 0;
    s = (1280 * rate) / 8000;
    return s + s / 8;
}

/* encoder.hh / decoder.hh fetch(): the polar code of each operation mode. */
const uint32_t *rb_mode_frozen(int mode, int *data_bits)
{
    switch (mode) {
    case 14:
        *data_bits = 1360;
        return rb_frozen_2048_1392;
    case 15:
        *data_bits = 1024;
        return rb_frozen_2048_1056;
    case 16:
        *data_bits = 680;
        return rb_frozen_2048_712;
    default:
        *data_bits = 0;
        return 0;
    }
}
