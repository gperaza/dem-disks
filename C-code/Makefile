GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CFLAGS = -Wall -O2 -std=gnu99 $(shell pkg-config --cflags  cairo)
CFLAGS += -DVERSION=\"$(GIT_VERSION)\"
CFLAGS += -DCOLLISIONS
#CFLAGS += -DGRAPHICS
GRAPHLIBS = $(shell pkg-config --libs cairo)
LIBS = -lm -lgsl -lgslcblas
CC = gcc

HDRS = common.h particle_structures.h
OBJS = main.o graphics.o integrator.o forces.o linkcell.o write.o

default:$(OBJS)
	$(CC)  $(CFLAGS) $(OBJS) $(GRAPHLIBS) $(LIBS) -o ../disks

clean:
	-/bin/rm	$(OBJS) disks


main.o:		Makefile $(HDRS)
integrator.o:	Makefile $(HDRS)
forces.o:	Makefile $(HDRS)
graphics.o:	Makefile $(HDRS)
linkcell.o:	Makefile $(HDRS)
write.o:	Makefile $(HDRS)
