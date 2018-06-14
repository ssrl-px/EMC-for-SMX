.PHONY: clean
CC = mpicc
CFLAGS = -O3
LFLAGS = -lm

SRC = recon.c setup.c tomo.c max.c Aval.c merge-sort.c
OBJ = $(patsubst %.c, %.o, $(SRC))
DEP = emc.h
EXE = emc

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LFLAGS) $(CFLAGS) -o $@

$(OBJ): %.o: %.c $(DEP)
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	rm *.o
