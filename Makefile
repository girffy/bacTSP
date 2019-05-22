WUSS=-Wno-unused-local-typedefs -Wno-unused-variable -Wno-unused-but-set-variable
CFLAGS+=-Wall -Werror -Wextra -Wpedantic ${WUSS}

debug:
	g++ -DDEBUG ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 #-pg -no-pie

release:
	g++ ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 -march=native #-pg -no-pie

clean:
	rm a.out

# # We need a test script
# make debug CFLAGS= && mv a.out debug_b2cut
# make debug CFLAGS=-DBOOSTCUT && mv a.out debug_boostcut
# make release CFLAGS= && mv a.out release_b2cut
# make release CFLAGS=-DBOOSTCUT && mv a.out release_boostcut
# for i in debug_b* release_b* ; do echo $i ; time ./$i 50 >/dev/null ; echo "" ; done
