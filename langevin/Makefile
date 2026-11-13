LINLIBS = GSL
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -g
DEPFLAGS := -MMD -MP
LDLIBS := -lm
LDFLAGS :=
ifeq ($(LINLIBS), LAPACK)
CXXFLAGS += -DLAPACK
LDLIBS += -lgsl -llapacke -lcblas -llapack -lblas -lgfortran
else ifeq ($(LINLIBS), MKL)
CXXFLAGS += -DMKL
LDLIBS += -lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#source /opt/intel/mkl/bin/mklvars.sh intel64
else
CXXFLAGS += -DGSL
LDLIBS += -lgsl -lgslcblas
endif
LDLIBS += -lyaml-cpp

SRC = main.cpp model.cpp utils.cpp forcefield.cpp neighborlist.cpp linalg.cpp stepper.cpp constraint.cpp integration.cpp yaml_config.cpp
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_langevin
PREFIX = $(CONDA_PREFIX)

all: $(PROG)

$(PROG): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

# tests: force/energy consistency and golden-run regression
TEST_OBJS = $(filter-out main.o, $(OBJS))
CHECK = ../bin/check_forces

$(CHECK): ../tests/check_forces.cpp $(TEST_OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -MF check_forces.d -MT $@ -I. $< $(TEST_OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY: all check golden clean install

check: $(CHECK)
	$(CHECK)

golden: $(PROG)
	bash ../tests/golden/compare.sh

clean:
	rm -f *~ *.o *.d $(PROG) $(CHECK) core a.out

install: $(PROG)
	cp $(PROG) $(PREFIX)/bin/

-include $(OBJS:.o=.d) check_forces.d
