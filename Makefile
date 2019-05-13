CFLAGS=-Wall -Werror -Wextra -Wpedantic -Wno-unused-local-typedefs

all:
	g++ ${CFLAGS} *.cpp -lglpk -std=c++1z -g -O3

clean:
	rm a.out
