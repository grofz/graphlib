# A dumb Makefile
# 2024
#
.SUFFIXES:

.PHONY: all clean directories

# path to sources and objects
DIR = build
BINDIR = bin
vpath %.f90 src:testsrc

# path to modules
JDIR = include

# compiler flags
COMMON_FLAGS = -Wall -Wextra -pedantic -std=gnu -fimplicit-none -fbacktrace -fPIC -fmax-errors=5 -cpp -MMD
# to provide debugging executable, run
#   $ make clean
#   $ make DEBUG=yes
#
ifdef DEBUG
  FFLAGS =-Og -g -fcheck=all -Wtrampolines $(COMMON_FLAGS)
else
  FFLAGS =-Ofast -march=native $(COMMON_FLAGS)
endif

# compiler and linker
FC = gfortran
AR = ar -rcv

# object files
# all .f90 files from the src/ directory (automatic list)
# 1. Find every .f90 file inside the src/ directory
SRC_FILES = $(wildcard src/*.f90)
# 2. Strip the 'src/' prefix and change the extension from '.f90' to '.o'
OBJ_NAMES = $(patsubst src/%.f90, %.o, $(SRC_FILES))
# 3. Add the 'build/' directory prefix to all of them
MODOBJECTS = $(addprefix $(DIR)/, $(OBJ_NAMES))

MAIN1 = build_test/vtuiotest.o
MAIN2 = build_test/maxflowtest.o
MAIN3 = build_test/betweennesstest.o
MAIN4 = build_test/concomtest.o
ALLOBJECTS = $(MODOBJECTS) $(MAIN1) $(MAIN2) $(MAIN3) $(MAIN4)

# dependent libraries
#ldir = ./lib/odepack
#libs = -lodepack

# output library
OUTLIB = libgraph.a

ifdef OS
	EXE1 = test_vtuio.exe
	EXE2 = test_maxflow.exe
	EXE3 = test_betweeness.exe
	EXE4 = test_concom.exe
else
	EXE1 = $(BINDIR)/test_vtuio.x
	EXE2 = $(BINDIR)/test_maxflow.x
	EXE3 = $(BINDIR)/test_betweeness.x
	EXE4 = $(BINDIR)/test_concom.x
endif

# default goal and dependencies
all: directories $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(OUTLIB)

# Ensure directories exist before compilation begins
directories:
	@mkdir -p $(DIR) build_test $(BINDIR) $(JDIR)

$(EXE1) : $(MODOBJECTS) $(MAIN1)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE2) : $(MODOBJECTS) $(MAIN2)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE3) : $(MODOBJECTS) $(MAIN3)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE4) : $(MODOBJECTS) $(MAIN4)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^

$(OUTLIB) : $(MODOBJECTS)
	$(AR) $@ $^

$(ALLOBJECTS) : Makefile

# module dependencies
$(DIR)/graph.o : $(DIR)/graph_user.o $(DIR)/graph_adjlist.o $(DIR)/conts.o

#$(DIR)/graph_smod_centrality.o : $(DIR)/graph.o
#$(DIR)/graph_smod_flow.o : $(DIR)/graph.o

$(DIR)/graph_testutils.o : $(DIR)/graph_user.o $(DIR)/graph.o $(DIR)/parse.o

$(DIR)/vtuio26.o : $(DIR)/graph.o $(DIR)/vtuio_tree.o

# dependencies for unit test programs
#
#build_test/vtuiotest.o : $(MODOBJECTS)

#build_test/maxflowtest.o : $(MODOBJECTS)

#build_test/betweenesstest.o : $(MODOBJECTS)

#build_test/concomtest.o : $(MODOBJECTS)

#.f.o:
build_test/%.o : testsrc/%.f90
	$(FC) $(FFLAGS) -J$(JDIR) -c $< -o $@
$(DIR)/%.o : src/%.f90
	$(FC) $(FFLAGS) -J$(JDIR) -c $< -o $@

# phony clean-up target
clean :
	-rm -f $(DIR)/*.o build_test/*.o $(JDIR)/*.mod $(JDIR)/*.smod $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(OUTLIB)

# Include the generated dependency files if they exist
-include $(MODOBJECTS:.o=.d) $(ALLOBJECTS:.o=.d)
# end of makefile
