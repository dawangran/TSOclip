CC ?= gcc
CFLAGS ?= -O3 -march=native -pipe -fopenmp -DNDEBUG
LDFLAGS ?= -lz

BIN = tsoclip
SRC = src/tsoclip.c

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(BIN)

.PHONY: all clean
