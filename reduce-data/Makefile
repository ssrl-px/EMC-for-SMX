.PHONY: clean
CC = gcc
CFLAGS = -O3
LFLAGS = -lm -lssl -lcrypto

SRC = write-data.c read-cbf.c
OBJ = $(patsubst %.c, %.o, $(SRC))
DEP = write-data.h
EXE = wr-data

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LFLAGS) $(CFLAGS) -o $@

$(OBJ): %.o: %.c $(DEP)
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	rm *.o
