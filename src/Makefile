IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -ggdb

ODIR=../obj
LDIR =../lib


LIBS=-lm

_DEPS = lbm.h utilities.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = lbm.o utilities.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

lbm: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *~ core $(INCDIR)/*~ 