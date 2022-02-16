CC = gcc
CFLAGS = -Wall -shared -fPIC

all: local bfx kalign
.PHONY: all

local:
	$(CC) $(CFLAGS) align.c -o libalign.so
	$(CC) $(CFLAGS) seqtools.c -o libseqtools.so
	$(CC) $(CFLAGS) consensus.c -o libconsensus.so
.PHONY: local

bfx:
	if [ -f bfx/Makefile ]; then make -C bfx; fi
.PHONY: bfx

kalign:
	if [ -f kalign/Makefile ]; then make -C kalign; fi
.PHONY: kalign

clean: clean_local clean_bfx clean_kalign
.PHONY: clean

clean_local:
	rm -f libalign.so libseqtools.so libconsensus.so bfx/libswalign.so
.PHONY: clean_local

clean_bfx:
	if [ -f bfx/Makefile ]; then make -C bfx clean; fi
.PHONY: clean_bfx

clean_kalign:
	if [ -f kalign/Makefile ]; then make -C kalign clean; fi
.PHONY: clean_kalign
