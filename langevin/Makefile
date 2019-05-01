CC = g++
LINLIBS = GSL
#CFLAGS = -std=c++11 -O3 -Wall
ifeq ($(LINLIBS), LAPACK)
CFLAGS = -std=c++11 -O3 -DLAPACK
#LDLIBS = -lgsl -lsbtk -llapacke -lcblas -llapack -lblas -lgfortran -lm
LDLIBS = -lgsl -llapacke -lcblas -llapack -lblas -lgfortran -lm
else ifeq ($(LINLIBS), MKL)
CFLAGS = -std=c++11 -O3 -DMKL
LDLIBS = -lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
#source /opt/intel/mkl/bin/mklvars.sh intel64
else
CFLAGS = -std=c++11 -O3 -DGSL
LDLIBS = -lgsl -lgslcblas -lm
endif

LDFLAGS =
SRC = main.cpp model.cpp utils.cpp forcefield.cpp linalg.cpp stepper.cpp
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_langevin

all: $(PROG)

$(PROG): $(OBJS)
	 $(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

main.o: main.cpp utils.h model.h
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
