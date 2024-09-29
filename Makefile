ARCH := $(shell $(CC) -dumpmachine | awk -F"-" '{print $$1}')

$(info $(ARCH))

TARGETS = generic generic_ossl

generic_SOURCES = keccak1600*.c sha3.c sha3_test.c
generic_CFLAGS := $(EXTRA_CFLAGS) -O2

generic_ossl_SOURCES = sha3_ossl.c sha3_test.c
generic_ossl_CFLAGS := $(EXTRA_CFLAGS) -O2 -DOSSL_BUILD
generic_ossl_LIBS = -lcrypto

ifeq ($(ARCH),x86_64)
	EXTRA_CFLAGS := -march=native -mtune=native
endif

ifeq ($(ARCH),riscv64)
	generic_SOURCES += keccak1600_intermediateur_rv64i.S keccak1600_inplaceur_rv64id.S
	generic_CFLAGS += -DRVASM_IMPL
endif

.PHONY: all clean

all: $(TARGETS) clean-objs

$(TARGETS):
	$(CC) -o sha3_$@ $($@_CFLAGS) $($@_SOURCES) $($@_LIBS)

clean-objs:
	rm -f *.o

clean: clean-objs
	rm -f $(foreach target, $(TARGETS), sha3_$(target))
