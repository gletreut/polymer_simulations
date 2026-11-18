CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -g
DEPFLAGS := -MMD -MP
LDLIBS := -lgsl -lgslcblas -lm
SRC = main.cpp model.cpp geometry.cpp utils.cpp
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_montecarlo
PREFIX = $(CONDA_PREFIX)

all: $(PROG)

$(PROG): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

.PHONY: all clean install

clean:
	rm -f *~ *.o *.d $(PROG) core a.out

install: $(PROG)
	cp $(PROG) $(PREFIX)/bin/

-include $(OBJS:.o=.d)
