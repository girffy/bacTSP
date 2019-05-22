CFLAGS=-Wall -Werror -Wextra -Wpedantic -Wno-unused-local-typedefs

all:
	g++ ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 #-pg -no-pie

clean:
	rm a.out
