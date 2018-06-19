.PHONY: clean
CC = mpicc
CFLAGS = -O3
LFLAGS = -lm -lssl -lcrypto

SRC = mpi-make-background.c read-cbf.c
OBJ = $(patsubst %.c, %.o, $(SRC))
DEP = make-bg.h
EXE = ave_bg

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LFLAGS) $(CFLAGS) -o $@

$(OBJ): %.o: %.c $(DEP)
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	rm *.o
