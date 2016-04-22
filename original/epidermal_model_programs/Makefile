#make
CC = icc 
CFLAGS = 
PROGRAM = epidermis

OBJS	= auto21-1.o 21-der.o cell_dynamics.o ca_dynamics.o func21-2.o

$(PROGRAM):	$(OBJS)
		$(CC) $(OBJS) -DMEXP=607 SFMT.c -lm -o $(PROGRAM)

auto21-1.o: main.h param.h funcs.h
21-der.o: param.h funcs.h
cell_dynamics.o: param.h funcs.h
ca_dynamics.o: param.h funcs.h
func21-2.o: param.h funcs.h

clean:
	rm -f *.o *~ \#*
