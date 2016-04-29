ARCH            := $(shell uname -s | tr a-z A-Z)

PROG	:=	sandpiles
BINDIR = ../bin
CC	= gcc
CFLAGS := -g  -O1 -std=c99 -Wno-deprecated-declarations -D DIM=$(DIM) 	\
-D CASE=$(CASE) -D METHOD=$(METHOD)

ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS		:= -lOpenCL -lGL -lGLU -lglut -lm -fopenmp -pthread
endif

.phony: default clean clear

default: $(PROG)

$(PROG): main.o display.o sand.o vector.o util.o
	mkdir -p $(BINDIR)
	$(CC) -o $(BINDIR)/$@-m$(METHOD)-c$(CASE)-d$(DIM)	\
$(LDFLAGS) $^ $(LDLIBS)

main.o: main.c display.h type.h sand.h vector.h util.h
	$(CC) $(CFLAGS) -c -o main.o main.c $(LDLIBS)

sand.o: sand.c sand.h type.h

display.o: display.c display.h type.h

vector.o: vector.c util.h type.h sand.h

util.o: util.h

clear:
	rm -f *.o

clean: clear
	rm -f $(BINDIR)/$(PROG)-* *~
