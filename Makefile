CFLAGS=-Wall -Werror -Wextra -Wpedantic -Wno-unused-local-typedefs

all:
	g++ ${CFLAGS} *.cpp -lglpk -std=c++2a

clean:
	rm a.out
