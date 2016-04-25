MAX_HEIGHT?= 4
DIM?= 128
CASE?= 1
METHOD?= 1

ARCH            := $(shell uname -s | tr a-z A-Z)

PROG	:=	sandpiles
BINDIR = ../bin
CC	= gcc
CFLAGS := -g -O3 -std=c99 -Wno-deprecated-declarations -D DIM=$(DIM)	\
-D MAX_HEIGHT=$(MAX_HEIGHT) -D CASE=$(CASE) -D METHOD=$(METHOD)
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL
else
LDLIBS		:= -lOpenCL -lGL -lGLU -lglut -lm -fopenmp
endif

.phony: default clean clear

default: $(PROG)

$(PROG): main.o display.o
	$(CC) -o $(BINDIR)/$@-m$(METHOD)-c$(CASE)-h$(MAX_HEIGHT)-d$(DIM) \
$(LDFLAGS) $^ $(LDLIBS)

main.o: main.c display.h
	$(CC) $(CFLAGS) -c -o main.o main.c $(LDLIBS)

display.o: display.c display.h

clear:
	rm -f *.o

clean: clear
	rm -f $(PROG)-* *~
