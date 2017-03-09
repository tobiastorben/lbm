CC=gcc
CFLAGS=-I$(IDIR) -std=c99 -m64 -Ofast -flto -march=native -funroll-loops

ODIR=../obj
LDIR =../lib
IDIR =../include

LIBS= -lm -lpthread

_DEPS = lbm.h utilities.h boundaries.h core.h initialization.h postProcessing.h inputParser.h structs.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = lbm.o utilities.o boundaries.o core.o initialization.o postProcessing.o inputParser.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

lbm: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	
clean:
	rm -f $(ODIR)/*
	rm lbm.exe