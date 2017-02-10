IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -std=c99 -m64 -Ofast -flto -march=native -funroll-loops

ODIR=../obj
LDIR =../lib


LIBS=

_DEPS = utilities.h boundaries.h core.h initialization.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = lbm.o utilities.o boundaries.o core.o initialization.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

lbm: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*