CC = g++
CFLAGS = -Wall -Werror -O3 -std=c++11
LIBS = -lm

LIBUQ = libuq.a
AR = ar
ARFLAGS = rv

SRCS = main.cpp helper.cpp
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

TARGET = UQ_twodim

all: $(TARGET)
	
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): main.cpp helper.cpp $(LIBUQ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBUQ) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o $(TARGET) $(LIBUQ)
