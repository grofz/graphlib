# A dumb Makefile
# 2024
#
.SUFFIXES:

.PHONY: all test clean directories

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
  FFLAGS =-Og -g -fcheck=all -Wtrampolines -DDEBUG $(COMMON_FLAGS)
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

TESTUTILS = build_test/map.o build_test/utils.o

MAIN1 = build_test/vtuiotest.o
MAIN2 = build_test/maxflowtest.o
MAIN3 = build_test/betweennesstest.o
MAIN4 = build_test/concomtest.o
MAIN5 = build_test/scctest.o
MAIN6 = build_test/contstest.o
MAIN_EXAMPLE = build_test/example.o

# output library
OUTLIB = libgraph.a

ifdef OS
	EXE1 = test_vtuio.exe
	EXE2 = test_maxflow.exe
	EXE3 = test_betweenness.exe
	EXE4 = test_concom.exe
	EXE5 = test_scc.exe
	EXE6 = test_conts.exe
	EXE_EXAMPLE = example.exe
else
	EXE1 = $(BINDIR)/test_vtuio
	EXE2 = $(BINDIR)/test_maxflow
	EXE3 = $(BINDIR)/test_betweenness
	EXE4 = $(BINDIR)/test_concom
	EXE5 = $(BINDIR)/test_scc
	EXE6 = $(BINDIR)/test_conts
	EXE_EXAMPLE = $(BINDIR)/example
endif

# default goal and dependencies
all: directories $(OUTLIB)
test: directories $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE_EXAMPLE) $(OUTLIB)

# Ensure directories exist before compilation begins
directories:
	@mkdir -p $(DIR) build_test $(BINDIR) $(JDIR)

$(EXE1) : $(MODOBJECTS) $(TESTUTILS) $(MAIN1)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE2) : $(MODOBJECTS) $(TESTUTILS) $(MAIN2)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE3) : $(MODOBJECTS) $(TESTUTILS) $(MAIN3)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE4) : $(MODOBJECTS) $(TESTUTILS) $(MAIN4)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE5) : $(MODOBJECTS) $(TESTUTILS) $(MAIN5)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE6) : $(MODOBJECTS) $(TESTUTILS) $(MAIN6)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^
$(EXE_EXAMPLE) : $(MODOBJECTS) $(MAIN_EXAMPLE)
	$(FC) $(FFLAGS) -J$(JDIR) -o $@ $^

$(OUTLIB) : $(MODOBJECTS)
	$(AR) $@ $^

# module dependencies (just for bootstraping)
build_test/vtuiotest.o : build_test/map.o

$(DIR)/graph.o : $(DIR)/graph_user.o $(DIR)/graph_adjlist.o $(DIR)/conts.o

build_test/utils.o : $(DIR)/graph_user.o $(DIR)/graph.o $(DIR)/parse.o

$(DIR)/vtuio26.o : $(DIR)/graph.o $(DIR)/vtuio_tree.o

$(DIR)/conts_test.o : $(DIR)/utest.o

#.f.o:
build_test/%.o : testsrc/%.f90
	$(FC) $(FFLAGS) -I$(JDIR) -Jbuild_test -c $< -o $@
$(DIR)/%.o : src/%.f90
	$(FC) $(FFLAGS) -J$(JDIR) -c $< -o $@

# phony clean-up target
clean :
	-rm -f $(DIR)/*.o build_test/*.o build_test/*.mod build_test/*.smod $(JDIR)/*.mod $(JDIR)/*.smod $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE_EXAMPLE) $(OUTLIB)

# Include the generated dependency files if they exist
-include $(MODOBJECTS:.o=.d) $(TESTUTILS:.o=.d)
# end of makefile
