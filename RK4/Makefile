CC=gcc
COPTS=-g -Wall -pedantic -Wextra -std=gnu99
ALL=RK4

all: $(ALL)

JUNK= *.o *~ *.dSYM RK4

clean: -rm -rf $(JUNK)

RK4: RK4.c
	$(CC) $(COPTS) $^ -lm -o $@
