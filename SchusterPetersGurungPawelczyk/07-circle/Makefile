#CC = scorep mpicc
CC = mpicc
CFLAGS = -std=c11 -pedantic -Wall -Wextra -Og -gdwarf-4
LFLAGS = $(CFLAGS)
TGTS = circle.x

all: $(TGTS)

$(TGTS): %.x: %.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(TGTS)
