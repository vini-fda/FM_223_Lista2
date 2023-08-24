CC = gcc
CFLAGS = -Wall -shared -fPIC -O3

LIB_NAME = libhenon.so
SRC_FILE = henon.c

all: $(LIB_NAME)

$(LIB_NAME): $(SRC_FILE)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f $(LIB_NAME)

.PHONY: clean
