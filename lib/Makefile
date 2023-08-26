CC = gcc
CFLAGS = -Wall -shared -fPIC -O3

all: libhenon.so libpendulum.so

libhenon.so: henon.c
	$(CC) $(CFLAGS) -o $@ $^

libpendulum.so: pendulum.c
	$(CC) $(CFLAGS) -lm -o $@ $^

clean:
	rm -f libhenon.so libpendulum.so

.PHONY: clean
