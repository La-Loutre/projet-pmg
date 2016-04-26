ARCH            := $(shell uname -s | tr a-z A-Z)

PROG	:=	sandpiles
BINDIR = ../bin
CC	= gcc
CFLAGS := -g -O1 -std=c99 -Wno-deprecated-declarations -D DIM=$(DIM) 	\
-D CASE=$(CASE) -D METHOD=$(METHOD)
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
	$(CC) -o $(BINDIR)/$@-m$(METHOD)-c$(CASE)-d$(DIM)	\
$(LDFLAGS) $^ $(LDLIBS)

main.o: main.c display.h
	$(CC) $(CFLAGS) -c -o main.o main.c $(LDLIBS)

display.o: display.c display.h

clear:
	rm -f *.o

clean: clear
	rm -f $(BINDIR)/$(PROG)-* *~
