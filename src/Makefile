#CC = mpicxx
CC = g++
CFLAGS = -Wall -Werror -O2
.cpp.o:  ; $(CC) -cpp $(CFLAGS) $<

OBJ = 	main.o\
		helper.o

all:  $(OBJ)
	$(CC) $(CFLAGS) -o UQ $(OBJ)  -lm

%.o : %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o

clean:
	rm $(OBJ) UQ

helper.o : helper.hpp

main.o : helper.hpp MCsim_normal.hpp MCsim_uniform.hpp
