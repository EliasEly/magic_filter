CC=clang
OPTIONS=-O3 -g -Wall -Werror -march=native
EXEC=magic_exec

all: $(EXEC)

magic_exec: magic_exec.o MagicFilter_naive.o
	$(CC) -o $@ $^

magic_exec.o: magic_exec.c magic_exec.h
	$(CC) $(OPTIONS) -o $@ -c $< 

MagicFilter_NEON.o: MagicFilter_NEON.c magic_exec.h
	$(CC) $(OPTIONS) -o $@ -c $<

MagicFilter_naive.o: MagicFilter_naive.c magic_exec.h
	$(CC) $(OPTIONS) -o $@ -c $<

MagicFilter_sse.o: MagicFilter_sse.c magic_exec.h
	$(CC) $(OPTIONS) -o $@ -c $<

clean:
	rm -f *.o $(EXEC)
