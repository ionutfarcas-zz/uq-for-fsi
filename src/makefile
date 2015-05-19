CC = g++
CFLAGS = -Wall -Werror -O2
LIBS = -lm

LIBUQ = libuq.a
AR = ar
ARFLAGS = rv

SRCS = main.cpp helper.cpp
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: UQ
	
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

UQ: main.cpp helper.cpp $(LIBUQ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBUQ) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o UQ $(LIBUQ)
