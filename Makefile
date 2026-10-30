# ribbit-core -- plain make, C99, no dependencies besides libm.
#   make            static + shared library into $(BUILD)
#   make test       build and run the C tests (goldens from the C++ reference, round trips)
#   make sanitize   the C tests under AddressSanitizer + UndefinedBehaviorSanitizer
#   make ref        build the C++ reference harness against ref/rattlegram and compare it
#                   with the core on fresh vectors (needs a C++17 compiler and python3)
#   make goldens    regenerate tests/golden/ with the C++ reference (needs ref/rattlegram)
#   make BUILD=/abs/dir shared   out-of-tree build (used by wave-owl)
BUILD ?= build
CC ?= cc
CXX ?= c++
CFLAGS ?= -O2
RATTLEGRAM ?= ref/rattlegram/app/src/main/cpp
# -ffp-contract=off: no fused multiply-add on either side. The reference is compiled with
# the same flag, so the float arithmetic follows the same expression order bit for bit.
CORE_FLAGS = -std=c99 -Wall -Wextra -Werror -pedantic -ffp-contract=off -fPIC -Iinclude
LDLIBS = -lm
TEST_FLAGS = -std=c99 -Wall -Wextra -Werror -pedantic -ffp-contract=off -Iinclude -Isrc
REF_FLAGS = -std=c++17 -O2 -ffp-contract=off -w -I$(RATTLEGRAM)

UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
  SHLIB = libribbit.dylib
  SHFLAGS = -dynamiclib -install_name @rpath/libribbit.dylib
else
  SHLIB = libribbit.so
  SHFLAGS = -shared -Wl,-soname,libribbit.so
endif

HDRS = include/ribbit.h src/ribbit_internal.h
OBJ = $(BUILD)/ribbit_core.o $(BUILD)/ribbit_fft.o $(BUILD)/ribbit_code.o \
      $(BUILD)/ribbit_tables.o $(BUILD)/ribbit_rx.o $(BUILD)/ribbit_tx.o

all: static shared
static: $(BUILD)/libribbit.a
shared: $(BUILD)/$(SHLIB)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.c $(HDRS) | $(BUILD)
	$(CC) $(CORE_FLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/libribbit.a: $(OBJ)
	ar rcs $@ $^

$(BUILD)/$(SHLIB): $(OBJ)
	$(CC) $(SHFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/test_ribbit: tests/test_ribbit.c $(BUILD)/libribbit.a
	$(CC) $(TEST_FLAGS) $(CFLAGS) $< $(BUILD)/libribbit.a -o $@ $(LDLIBS)

$(BUILD)/ribbit_tool: tools/ribbit_tool.c $(BUILD)/libribbit.a
	$(CC) $(TEST_FLAGS) $(CFLAGS) $< $(BUILD)/libribbit.a -o $@ $(LDLIBS)

$(BUILD)/ref_tool: tools/ref_tool.cpp | $(BUILD)
	@test -f $(RATTLEGRAM)/decoder.hh || { echo "ref/rattlegram missing: git submodule update --init ref/rattlegram"; exit 1; }
	$(CXX) $(REF_FLAGS) $< -o $@

test: $(BUILD)/test_ribbit
	$(BUILD)/test_ribbit tests/golden

SAN_FLAGS = -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all

sanitize:
	$(MAKE) --no-print-directory BUILD=$(BUILD)/san CFLAGS="$(SAN_FLAGS)" \
		LDLIBS="$(LDLIBS) $(SAN_FLAGS)" test

tool: $(BUILD)/ribbit_tool

ref: $(BUILD)/ribbit_tool $(BUILD)/ref_tool
	python3 tools/ref_check.py --core $(BUILD)/ribbit_tool --ref $(BUILD)/ref_tool --work $(BUILD)/refcheck

goldens: $(BUILD)/ref_tool
	python3 tools/ref_check.py --ref $(BUILD)/ref_tool --goldens tests/golden

clean:
	rm -rf $(BUILD)

.PHONY: all static shared test sanitize tool ref goldens clean
