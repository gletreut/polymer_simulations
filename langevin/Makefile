CC = g++
#CFLAGS = -std=c++11 -O3 -Wall
CFLAGS = -std=c++11 -O3
LDLIBS = -lgsl -lm
LDFLAGS =
SRC = main.cpp model.cpp utils.cpp
OBJS = $(SRC:.cpp=.o)
PROG = ../bin/simu_langevin

all: $(PROG)

$(PROG): $(OBJS)
	 $(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

main.o: main.cpp utils.h model.h
	$(CC) $(CFLAGS) -c main.cpp

model.o: model.h model.cpp
	$(CC) $(CFLAGS) -c model.cpp

utils.o: utils.h utils.tpp utils.cpp
	$(CC) $(CFLAGS) -c utils.cpp

clean:
	rm -f *~ *.o $(PROG) core a.out
