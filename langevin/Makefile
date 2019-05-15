CC = g++
LINLIBS = GSL
CFLAGS := -std=c++11 -O3 -Wall
LDLIBS := -lm
LDFLAGS :=
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

SRC = main.cpp model.cpp utils.cpp forcefield.cpp linalg.cpp stepper.cpp
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_langevin

all: $(PROG)

$(PROG): $(OBJS)
	 $(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

main.o: main.cpp utils.h linalg.h forcefield.h model.h stepper.h global.h
	$(CC) $(CFLAGS) -c main.cpp

model.o: model.h model.cpp
	$(CC) $(CFLAGS) -c model.cpp

forcefield.o: forcefield.h forcefield.cpp
	$(CC) $(CFLAGS) -c forcefield.cpp

stepper.o: stepper.h stepper.cpp
	$(CC) $(CFLAGS) -c stepper.cpp

linalg.o: linalg.h linalg.cpp
	$(CC) $(CFLAGS) -c linalg.cpp

utils.o: utils.h utils.tpp utils.cpp
	$(CC) $(CFLAGS) -c utils.cpp

clean:
	rm -f *~ *.o $(PROG) core a.out
