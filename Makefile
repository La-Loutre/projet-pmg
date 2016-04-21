DIM ?= 128
MAX_HEIGHT ?= 4

ARCH            := $(shell uname -s | tr a-z A-Z)

PROG	:=	sable

CC	= gcc
CFLAGS := -g -O3 -std=c99 -Wno-deprecated-declarations -D DIM=$(DIM)	\
-D MAX_HEIGHT=$(MAX_HEIGHT)
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS		:= -lOpenCL -lGL -lGLU -lglut -lm
endif

.phony: default clean

default: case1 case2

case1: main_case1.o display.o
	$(CC) -o $(PROG)-$@-$(DIM)-$(MAX_HEIGHT) $(LDFLAGS) $^ $(LDLIBS)

case2: main_case2.o display.o
	$(CC) -o $(PROG)-$@-$(DIM)-$(MAX_HEIGHT) $(LDFLAGS) $^ $(LDLIBS)

main_case1.o: display.h
	$(CC) -D CASE=1 $(CFLAGS) -c -o main_case1.o main.c

main_case2.o: display.h
	$(CC) -D CASE=2 $(CFLAGS) -c -o main_case2.o main.c

#main.o: display.h

display.o: display.h

clean:
	rm -rf *.o $(PROG)-*
