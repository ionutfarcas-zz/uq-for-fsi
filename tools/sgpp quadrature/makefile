CC = g++
INC = -I ~/Tools/sgpp_trunk/base/src
CFLAGS = -O3 -std=c++11 $(INC)
LIBS = -lm -L ~/Tools/sgpp_trunk/lib/sgpp -lsgppbase


TARGET = test_sgpp

SRCS = sgpp_quadrature.cpp
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: $(TARGET)

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(SRCS) $(LIBOPT)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -rf *.o $(TARGET)

