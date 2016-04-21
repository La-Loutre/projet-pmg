
ARCH            := $(shell uname -s | tr a-z A-Z)

PROG	:=	sable

CFLAGS	:=	-g -O3 -std=c99 -Wno-deprecated-declarations
ifeq ($(ARCH),DARWIN)
CFLAGS	+=	-I /opt/local/include
LDFLAGS	+=	-L /opt/local/include
LDLIBS	+=	-framework GLUT -framework OpenGL -framework OpenCL
else
LDLIBS		:= -lOpenCL -lGL -lGLU -lglut -lm
endif

.phony: default clean

default: $(PROG)

$(PROG): main.o display.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

main.o: display.h

display.o: display.h

clean:
	rm -rf *.o $(PROG)
