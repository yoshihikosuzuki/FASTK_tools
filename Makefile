CC = gcc
CFLAGS = -g -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDFLAGS = -lpthread -lm
ALL = CoverRate

all: $(ALL)

libfastk.c: gene_core.c
libfastk.h: gene_core.h

CoverRate: CoverRate.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o $@ CoverRate.c libfastk.c $(LDFLAGS)

Homex: Homex.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o $@ Homex.c libfastk.c $(LDFLAGS)

HomexDi: HomexDi.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o $@ HomexDi.c libfastk.c $(LDFLAGS)

HomexTri: HomexTri.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o $@ HomexTri.c libfastk.c $(LDFLAGS)

DropContext: DropContext.c libfastk.c libfastk.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o $@ DropContext.c libfastk.c DB.c QV.c $(LDFLAGS)

clean:
	rm -f $(ALL)
