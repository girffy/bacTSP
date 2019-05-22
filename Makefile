WUSS=-Wno-unused-local-typedefs -Wno-unused-variable -Wno-unused-but-set-variable
CFLAGS=-Wall -Werror -Wextra -Wpedantic ${WUSS}

debug:
	g++ -DDEBUG ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 #-pg -no-pie

release:
	g++ ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 #-pg -no-pie

clean:
	rm a.out
