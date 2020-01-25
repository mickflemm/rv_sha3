ARCH := $(shell $(CC) -dumpmachine | awk -F"-" '{print $$1}')

$(info $(ARCH))

TARGETS = generic generic_ossl

ifeq ($(ARCH),x86_64)
	EXTRA_CFLAGS := -march=native -mtune=native
endif

ifeq ($(ARCH),riscv64)
	generic_SOURCES = keccak1600_rv64i.S
	TARGETS += rv64i
endif

generic_SOURCES = keccak1600.c keccak1600_sponge.c sha3.c sha3_test.c
generic_CFLAGS := $(EXTRA_CFLAGS) -O3 -DUSE_UNROLLED_CF

generic_ossl_SOURCES = sha3_ossl.c sha3_test.c
generic_ossl_CFLAGS := $(EXTRA_CFLAGS) -O3
generic_ossl_LIBS = -lcrypto

rv64i_SOURCES = keccak1600_rv64i.S keccak1600_sponge.c sha3.c sha3_test.c
rv64i_CFLAGS := -march=rv64i -mabi=lp64 -O3

.PHONY: all clean

all: $(TARGETS) clean-objs

$(TARGETS):
	$(CC) -o sha3_$@ $($@_CFLAGS) $($@_SOURCES) $($@_LIBS)

clean-objs:
	rm -f *.o

clean: clean-objs
	rm -f $(foreach target, $(TARGETS), sha3_$(target))
