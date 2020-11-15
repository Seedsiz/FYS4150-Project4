all: compile execute

compile:
	Clang++ -Wall -o main.out $(wildcard *.cpp) -larmadillo -lomp -std=c++11 -Ofast

execute:
	./main.out
