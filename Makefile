CC ?= cc
CPPFLAGS ?=
CFLAGS ?= -O3 -march=native -pipe -DNDEBUG
LDFLAGS ?= -lz
OPENMP ?= auto
PREFIX ?= /usr/local
VERSION := $(shell cat VERSION 2>/dev/null || printf "0.2.0-dev")
CPPFLAGS += -DTSOCLIP_VERSION=\"$(VERSION)\"

ifeq ($(OPENMP),auto)
OPENMP_OK := $(shell $(CC) -x c -fopenmp -c /dev/null -o /dev/null >/dev/null 2>&1 && echo 1)
else ifeq ($(OPENMP),1)
OPENMP_OK := 1
else
OPENMP_OK :=
endif

ifeq ($(OPENMP_OK),1)
CFLAGS += -fopenmp
endif

BIN = tsoclip
SRC = src/tsoclip.c

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(BIN)

check: $(BIN)
	sh tests/run_tests.sh ./$(BIN)

smoke: check

install: $(BIN)
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 "$(BIN)" "$(DESTDIR)$(PREFIX)/bin/$(BIN)"

uninstall:
	rm -f "$(DESTDIR)$(PREFIX)/bin/$(BIN)"

.PHONY: all clean check smoke install uninstall
