CC = gcc
XX = g++
AR = ar

CFLAGS = -O3 -Wall
#CPPFLAGS = -O3 -Wall -std=c++11
CPPFLAGS = -g -Wall -std=c++11

#C_FILE = $(wildcard *.c)
#CPP_FILE = $(wildcard *.cpp)

C_FILE = $(shell echo ./src/laspack/*.c)
CPP_FILE = $(shell echo ./src/*.cpp)


C_OBJ= $(C_FILE:%.c=%.o)
CPP_OBJ= $(CPP_FILE:%.cpp=%.o)

LASPACK_LIB = ./lib/

laspack: $(C_OBJ)
	$(AR) -r liblaspack.a $(C_OBJ) 
	mv ./liblaspack.a ./lib/
	

clean:
	rm -f ./src/laspack/*.o
	rm -f ./src/*.o


cycas2:  $(CPP_OBJ)
	$(XX) -o $@ $(CPP_OBJ) -L$(LASPACK_LIB) -llaspack $(CPPFLAGS)

.cpp.o:
	$(XX) -c $(CPPFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<


