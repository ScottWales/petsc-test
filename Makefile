# Makefile for Fortran automatic dependency generation
# ====================================================
#
# Author: Scott Wales
# 
# Copyright 2014 ARC Centre of Excellence for Climate Systems Science
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

all: check bin/plex

.PHONY:all clean check
.SUFFIXES:

# Tested with gfortran-4.8 and ifort
FC         = mpif90
CC         = mpicc
LD         = ${CC}

CPPFLAGS  += -Iinclude

# These locations are defaults for vagrant/travis
PFUNIT    ?= ${HOME}/pfunit
PETSC_DIR ?= ${HOME}/petsc

FCFLAGS   += -I$(PFUNIT)/mod
VPATH     += $(PFUNIT)/mod
PFPARSE    = $(PFUNIT)/bin/pFUnitParser.py

CFLAGS    += -I${PETSC_DIR}/include
FCFLAGS   += -I${PETSC_DIR}/include
VPATH     += ${PETSC_DIR}/include
LDFLAGS   += -L${PETSC_DIR}/lib
LDLIBS    += -lpetsc

CFLAGS    += -std=c99 -Wall -Werror

# Compiler detection
ifeq ($(findstring gcc,$(shell $(FC) -v 2>&1)),gcc)
    COMPILER_TYPE =  gnu

    FCFLAGS      += -g -fbacktrace
    FCFLAGS      += -fimplicit-none
    FCFLAGS      += -Wall -Wextra -Werror
    FCFLAGS      += -Iinclude -Jmod

    TESTFCFLAGS  += -Wno-unused
    TESTFCFLAGS  += -Wno-uninitialized
    TESTFCFLAGS  += -Wno-unused-parameter

else ifeq ($(findstring ifort,$(shell $(FC) -v 2>&1)),ifort)
    COMPILER_TYPE =  intel

CFLAGS += -check-pointers=rw
LDFLAGS += -check-pointers=rw
    FCFLAGS      += -g -traceback
    FCFLAGS      += -warn all -warn errors -check all
    FCFLAGS      += -Iinclude -module mod

    TESTFCFLAGS  += -Wno-unused-parameter
endif

# .mod files are stored in this directory
VPATH   += mod

# Get source files
SRC     := $(shell find src -name '*.f90' -type f)
TESTSRC := $(shell find src -name '*.pf' -type f)

# Get list of tests to run
TESTS    = $(patsubst src/%.pf,test/%,$(TESTSRC))
TESTFLAGS = -no-signal-handler

# Read dependencies
DEPS += $(patsubst %.f90,deps/%.d,$(SRC))
DEPS += $(patsubst src/%.pf,deps/obj/%.d,$(TESTSRC))
-include $(DEPS)

# Compile programs found by the dependency generation
all: $(FPROGRAMS) bin/plex

bin/benchmark2: obj/readmesh.o

# Run all tests
check: $(TESTS)
	@for test in $^; do echo -e "\n$$test"; mpirun -n 2 ./$$test ${TESTFLAGS}; done

# Cleanup
clean:
	$(RM) -r bin test obj deps mod

# Compile source files
obj/%.o: src/%.f90
	@mkdir -p $(dir $@)
	@mkdir -p mod
	$(FC) $(FCFLAGS) -c -o $@ $<
obj/%.o: src/%.c
	@mkdir -p $(dir $@)
	${CC} ${CPPFLAGS} ${CFLAGS} -c -o $@ $<


# Process pFunit tests
obj/%.F90: src/%.pf
	@mkdir -p $(dir $@)
	$(PFPARSE) $< $@

# Compile tests
obj/%.o: obj/%.F90
	@mkdir -p $(dir $@)
	@mkdir -p mod
	$(FC) $(FCFLAGS) ${TESTFCFLAGS} -c -o $@ $<

# Secondexpansion to calculate prerequisite modules
.SECONDEXPANSION:

# Link programs (Fortran programs should use ${FC})
${FPROGRAMS}: LD=${FC}
bin/%: obj/%.o $$(objectrequirements_%.o)
	@mkdir -p $(dir $@)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Link tests with driver
$(TESTS):$(PFUNIT)/include/driver.F90
test/%: obj/%.o $$(objectrequirements_obj/%.o)
	@mkdir -p $(dir $@)
	$(FC) $(FCFLAGS) $(TESTFCFLAGS) $(LDFLAGS) -L$(PFUNIT)/lib -DUSE_MPI -openmp -DSUITE=$(or $(basename $(modulesprovided_obj/$*.o)),$(notdir $*))_suite -o $@ $^ $(LDLIBS) -lpfunit

# Dependency generation
deps/%.d: %.f90
	@mkdir -p $(dir $@)
	./gendeps -o $@ $<
deps/%.d: %.F90
	@mkdir -p $(dir $@)
	./gendeps -o $@ $<

