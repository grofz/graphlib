# A dumb Makefile
# 2024
#
.SUFFIXES:

.PHONY: all clean

# path to sources and objects
DIR = build
BINDIR = bin
vpath %.f90 src:test
vpath %.o $(DIR) build_test

# path to modules
JDIR = include

# compiler flags
# to provide debugging executable, run
#   $ make clean
#   $ make DEBUG=yes
#
ifdef DEBUG
  FFLAGS =-Og -g -Wall -Wextra -pedantic -std=gnu -fimplicit-none -fcheck=all -fbacktrace -fPIC
else
  FFLAGS =-Ofast -Wall -Wextra -pedantic -std=gnu -fimplicit-none -fbacktrace -fPIC
endif

# compiler and linker
FC = gfortran
AR = ar -rcv

# object files
MODOBJECTS = \
						 $(DIR)/conts.o \
						 $(DIR)/graph_adjlist.o \
						 $(DIR)/graph.o \
						 $(DIR)/vtuio_tree.o \
						 $(DIR)/vtuio26.o

MAIN1 = build_test/vtuiotest.o
MAIN2 = build_test/maxflowtest.o
ALLOBJECTS = $(MODOBJECTS) $(MAIN1) $(MAIN2)

# dependent libraries
#ldir = ./lib/odepack
#libs = -lodepack

# output library
OUTLIB = libgraph.a

ifdef OS
	EXE1 = test_vtuio.exe
	EXE2 = test_maxflow.exe
else
	EXE1 = $(BINDIR)/test_vtuio.x
	EXE2 = $(BINDIR)/test_maxflow.x
endif

# default goal and dependencies
all: $(EXE1) $(EXE2) $(OUTLIB)

$(EXE1) : $(MODOBJECTS) $(MAIN1)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE2) : $(MODOBJECTS) $(MAIN2)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^

$(OUTLIB) : $(MODOBJECTS)
	$(AR) $@ $^

$(ALLOBJECTS) : Makefile

# module dependencies
$(DIR)/graph.o : graph_adjlist.o conts.o

$(DIR)/vtuio26.o : graph.o vtuio_tree.o

# dependencies for unit test programs
#
build_test/vtuiotest.o : $(MODOBJECTS)

build_test/maxflowtest.o : $(MODOBJECTS)

#.f.o:
build_test/%.o : test/%.f90
	$(FC) $(FFLAGS) -J$(JDIR) -c $< -o $@
$(DIR)/%.o : src/%.f90
	$(FC) $(FFLAGS) -J$(JDIR) -c $< -o $@

# one phoney target
clean :
	-rm -f $(DIR)/*.o build_test/*.o $(JDIR)/*.mod $(EXE1) $(EXE2) $(OUTLIB)
# end of makefile
