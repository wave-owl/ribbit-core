/*
 * ref_tool -- the C++ reference: Rattlegram's own Encoder<RATE> and Decoder<RATE>
 * (ref/rattlegram/app/src/main/cpp, unmodified) behind the command line interface of
 * tools/ribbit_tool.c. The decoder is driven like Rattlegram's app (MainActivity.java):
 * rate/50 samples per feed(), process() when feed() returns true, staged() on SYNC, PING
 * and NOPE, fetch() on DONE.
 *
 *   ref_tool tx RATE CARRIER NOISE FANCY CALL PAYLOAD_HEX OUT.raw   (PAYLOAD_HEX "-": ping)
 *   ref_tool rx RATE IN.raw
 *
 * Copyright (c) 2026 Gerolf Ziegenhain. 0BSD, see include/ribbit.h.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <new>
#include <algorithm>
#define assert(expr) do {} while (0)
#include "encoder.hh"
#include "decoder.hh"

static int unhex(const char *hex, uint8_t *out, int max)
{
	int n = 0;
	if (!strcmp(hex, "-"))
		return 0;
	while (hex[0] && hex[1] && n < max) {
		unsigned v;
		if (sscanf(hex, "%2x", &v) != 1)
			return -1;
		out[n++] = (uint8_t)v;
		hex += 2;
	}
	return hex[0] ? -1 : n;
}

static EncoderInterface *make_encoder(int rate)
{
	switch (rate) {
	case 8000: return new Encoder<8000>();
	case 16000: return new Encoder<16000>();
	case 32000: return new Encoder<32000>();
	case 44100: return new Encoder<44100>();
	case 48000: return new Encoder<48000>();
	}
	return nullptr;
}

static DecoderInterface *make_decoder(int rate)
{
	switch (rate) {
	case 8000: return new Decoder<8000>();
	case 16000: return new Decoder<16000>();
	case 32000: return new Decoder<32000>();
	case 44100: return new Decoder<44100>();
	case 48000: return new Decoder<48000>();
	}
	return nullptr;
}

static int symbol_samples(int rate)
{
	int s = (1280 * rate) / 8000;
	return s + s / 8;
}

static int do_tx(int argc, char **argv)
{
	if (argc != 9)
		return 2;
	int rate = atoi(argv[2]);
	uint8_t payload[256] = {0};
	int len = unhex(argv[7], payload, 170);
	EncoderInterface *enc = make_encoder(rate);
	if (!enc || len < 0)
		return 1;
	int8_t call[16] = {0};
	strncpy(reinterpret_cast<char *>(call), argv[6], 9);
	enc->configure(payload, call, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]) != 0);
	std::vector<int16_t> buf(symbol_samples(rate));
	FILE *f = fopen(argv[8], "wb");
	if (!f)
		return 1;
	while (enc->produce(buf.data(), 0))
		fwrite(buf.data(), sizeof(int16_t), buf.size(), f);
	fclose(f);
	delete enc;
	return 0;
}

static void print_event(const char *kind, int mode, const char *call, float cfo, int flips,
			const uint8_t *payload, int len, long long sample)
{
	printf("%s mode=%d call=[%s] cfo=%a flips=%d len=%d sample=%lld payload=", kind, mode,
	       call, (double)cfo, flips, len, sample);
	for (int i = 0; i < len; ++i)
		printf("%02x", payload[i]);
	printf("\n");
}

static int do_rx(int argc, char **argv)
{
	if (argc != 4)
		return 2;
	int rate = atoi(argv[2]);
	DecoderInterface *dec = make_decoder(rate);
	FILE *f = fopen(argv[3], "rb");
	if (!dec || !f)
		return 1;
	std::vector<int16_t> all;
	int16_t tmp[4096];
	size_t n;
	while ((n = fread(tmp, sizeof(int16_t), 4096, f)) > 0)
		all.insert(all.end(), tmp, tmp + n);
	fclose(f);
	const int chunk = rate / 50;
	float sync_cfo = 0;
	char sync_call[10] = {0};
	int sync_mode = 0;
	for (size_t pos = 0; pos + chunk <= all.size(); pos += chunk) {
		long long consumed = (long long)(pos + chunk);
		if (!dec->feed(all.data() + pos, chunk, 0))
			continue;
		int status = dec->process();
		float cfo;
		int32_t mode;
		uint8_t call[10] = {0};
		uint8_t payload[171] = {0};
		switch (status) {
		case STATUS_SYNC:
		case STATUS_PING:
		case STATUS_NOPE:
			dec->staged(&cfo, &mode, call);
			print_event(status == STATUS_SYNC ? "SYNC" : status == STATUS_PING ? "PING" : "NOPE",
				    mode, reinterpret_cast<char *>(call), cfo, 0, payload, 0, consumed);
			if (status == STATUS_SYNC) {
				sync_cfo = cfo;
				sync_mode = mode;
				memcpy(sync_call, call, 10);
			}
			break;
		case STATUS_DONE: {
			int result = dec->fetch(payload);
			if (result < 0) {
				print_event("FAIL", sync_mode, sync_call, sync_cfo, -1, payload, 0, consumed);
			} else {
				int len = 0;
				while (len < 170 && payload[len])
					++len;
				print_event("DONE", sync_mode, sync_call, sync_cfo, result, payload, len, consumed);
			}
			break;
		}
		default:
			break;
		}
	}
	delete dec;
	return 0;
}

int main(int argc, char **argv)
{
	if (argc >= 2 && !strcmp(argv[1], "tx"))
		return do_tx(argc, argv);
	if (argc >= 2 && !strcmp(argv[1], "rx"))
		return do_rx(argc, argv);
	fprintf(stderr, "usage: ref_tool tx|rx ...\n");
	return 2;
}
