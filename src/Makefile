CC=gcc
CFLAGS=-I$(IDIR) -std=c99 -ggdb # -m64 -Ofast -flto -march=native -funroll-loops

ODIR=../obj
LDIR =../lib
IDIR =../include

LIBS= -lm

_DEPS = lbm.h boundaries.h core.h initialization.h inputParser.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = lbm.o boundaries.o core.o initialization.o inputParser.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

lbm: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	
clean:
	rm -f $(ODIR)/*