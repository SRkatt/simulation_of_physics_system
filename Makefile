FLAGS=-c -O2
LIB=-lm
NROBJ=nrutil.o jacobi.o lubksb.o ludcmp.o

all: 8.ps

8: 8.o $(NROBJ)
	gcc $(LIB) 8.o $(NROBJ) -o $@

8.out: 8
	./8

8.ps: 8.out 8.plt
	gnuplot 8.plt

.c.o:
	gcc $(FLAGS) $< -o $@

clean:
	rm -f *.o

distclean: clean
	rm -f 8 8.out 8.ps

