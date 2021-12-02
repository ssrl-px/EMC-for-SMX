.PHONY: clean
CC = mpicc
CFLAGS = -O3
LFLAGS = -lm -lssl -lcrypto -lcbf
INC = -I/home/sw/deb/CBFlib-0.9.6/include/

SRC = mpi-make-background.c
OBJ = $(patsubst %.c, %.o, $(SRC))
DEP = make-bg.h
EXE = ave_bg

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(INC)  $^ $(LFLAGS) $(CFLAGS)  -o $@

$(OBJ): %.o: %.c $(DEP)
	$(CC) $(INC)  -c $< $(CFLAGS) -o $@

clean:
	rm *.o
