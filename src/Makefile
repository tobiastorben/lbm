CFLAGS= -I$(IDIR) -std=c99 -m64 -Ofast -flto -march=native -funroll-loops -Wall -Wextra

ODIR=../obj
IDIR =../include

LIBS= -lm -lpthread

_DEPS = lbm.h output.h boundaries.h core.h initialization.h inputParser.h structs.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = lbm.o output.o boundaries.o core.o initialization.o inputParser.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	gcc -c -o $@ $< $(CFLAGS)

lbm: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
	
clean:
	rm -f $(ODIR)/*
