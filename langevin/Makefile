CC = g++
# CC = /bin/g++
LINLIBS = GSL
CFLAGS := -std=c++11 -O3 -Wall -g
LDLIBS := -lm
#LDFLAGS :=
LDFLAGS := -Wl,-rpath,$(CONDA_PREFIX)/lib/
ifeq ($(LINLIBS), LAPACK)
CFLAGS += -DLAPACK
#LDLIBS = -lgsl -lsbtk -llapacke -lcblas -llapack -lblas -lgfortran -lm
LDLIBS += -lgsl -llapacke -lcblas -llapack -lblas -lgfortran
else ifeq ($(LINLIBS), MKL)
CFLAGS += -DMKL
LDLIBS += -lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#source /opt/intel/mkl/bin/mklvars.sh intel64
else
CFLAGS += -DGSL
LDLIBS += -lgsl -lgslcblas
endif
LDLIBS += -lboost_filesystem -lboost_system
LDLIBS += -lyaml-cpp

SRC = main.cpp model.cpp utils.cpp forcefield.cpp neighborlist.cpp linalg.cpp stepper.cpp constraint.cpp integration.cpp yaml_config.cpp
# SRC = $(wildecard *.cpp)
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_langevin
PREFIX=$(CONDA_PREFIX)

all: $(PROG)

$(PROG): $(OBJS)
	 $(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

main.o: main.cpp utils.h linalg.h forcefield.h neighborlist.h model.h stepper.h global.h integration.h
	$(CC) $(CFLAGS) -c main.cpp

model.o: model.h forcefield.h utils.h constraint.h neighborlist.h model.cpp
	$(CC) $(CFLAGS) -c model.cpp

forcefield.o: forcefield.h linalg.h utils.h forcefield.cpp
	$(CC) $(CFLAGS) -c forcefield.cpp

constraint.o: constraint.h linalg.h utils.h constraint.cpp
	$(CC) $(CFLAGS) -c constraint.cpp

neighborlist.o: neighborlist.h linalg.h utils.h neighborlist.cpp
	$(CC) $(CFLAGS) -c neighborlist.cpp

integration.o: integration.h integration.cpp
	$(CC) $(CFLAGS) -c integration.cpp

stepper.o: stepper.h model.h linalg.h utils.h stepper.cpp
	$(CC) $(CFLAGS) -c stepper.cpp

linalg.o: linalg.h linalg.cpp
	$(CC) $(CFLAGS) -c linalg.cpp

utils.o: utils.h utils.tpp utils.cpp
	$(CC) $(CFLAGS) -c utils.cpp

yaml_config.o: yaml_config.h model.h stepper.h integration.h forcefield.h constraint.h neighborlist.h yaml_config.cpp
	$(CC) $(CFLAGS) -c yaml_config.cpp

.PHONY: clean
clean:
	rm -f *~ *.o $(PROG) core a.out

.PHONY: install
install: $(PROG)
	  cp $(PROG) $(PREFIX)/bin/
