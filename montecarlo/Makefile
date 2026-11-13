CXX ?= g++
CXXFLAGS = -std=c++11 -O3 -Wall
LDLIBS = -lgsl -lgslcblas -lm
LDFLAGS =
SRC = main.cpp model.cpp geometry.cpp utils.cpp
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_montecarlo

all: $(PROG)

$(PROG): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

main.o: main.cpp model.h utils.h
model.o: model.cpp model.h geometry.h
geometry.o: geometry.cpp geometry.h
utils.o: utils.cpp utils.h

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *~ *.o $(PROG) core a.out
