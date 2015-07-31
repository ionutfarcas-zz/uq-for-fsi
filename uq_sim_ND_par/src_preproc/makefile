CC = g++
INC = -I ~/Tools/sgpp_trunk/base/src
CFLAGS = -O3 -Werror -std=c++11 $(INC)
LIBS = -lm -L ~/Tools/sgpp_trunk/lib/sgpp -lsgppbase

LIBUQ = libuq.a
AR = ar
ARFLAGS = rv

SRCS = main.cpp helper.cpp
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

TARGET = UQ_ndim_par

all: $(TARGET)
	
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): main.cpp helper.cpp $(LIBUQ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBUQ) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o $(TARGET) $(LIBUQ)
