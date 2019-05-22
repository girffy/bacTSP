WUSS=-Wno-unused-local-typedefs -Wno-unused-variable -Wno-unused-but-set-variable
CFLAGS+=-Wall -Werror -Wextra -Wpedantic ${WUSS}

debug:
	g++ -DDEBUG ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 #-pg -no-pie

release:
	g++ ${CFLAGS} Graph.cpp TSP_prob.cpp main.cpp -lglpk -std=c++1z -g -O3 -march=native #-pg -no-pie

clean:
	rm a.out

# # We need a test script
bench:
	# just believe me when i say printf makes things rly rly slow
	# make debug   CFLAGS="-o debug_b2cut"
	# make debug   CFLAGS="-o debug_boostcut -DBOOSTCUT"
	make release CFLAGS="-o release_b2cut"
	make release CFLAGS="-o release_boostcut -DBOOSTCUT"

runbench:
	for i in release_b* ; do echo $$i ; time ./$$i 50 >/dev/null ; echo ; done

cleanbench:
	rm elease_b2cut release_boostcut
