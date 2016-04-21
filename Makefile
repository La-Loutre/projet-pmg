DIM?= 128
MAX_HEIGHT?= 4
CASE?= 1

ARCH := $(shell uname -s | tr a-z A-Z)
PROG :=	sandpiles
CC = gcc
CFLAGS := -g -O3 -std=c99 -Wno-deprecated-declarations -D DIM=$(DIM)	\
-D MAX_HEIGHT=$(MAX_HEIGHT) -D CASE=$(CASE)

ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS		:= -lOpenCL -lGL -lGLU -lglut -lm -fopenmp
endif

.PHONY: default clean clear

default:

seq: main.o display.o
	$(CC) -o $(PROG)-$@-case$(CASE)-$(DIM)-$(MAX_HEIGHT) $(LDFLAGS) \
$^ $(LDLIBS)

main.o: display.h

display.o: display.h

clear:
	rm -f *.o *~

clean: clear
	rm -f $(PROG)-*
