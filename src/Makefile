IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -ggdb -std=c99

ODIR=../obj
LDIR =../lib


LIBS=-lm

_DEPS = utilities.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = lbm.o utilities.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

lbm: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*