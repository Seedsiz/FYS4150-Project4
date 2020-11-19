all: compile execute

compile:
	g++-10 -Wall -o main.out $(wildcard *.cpp) -larmadillo -fopenmp -std=c++11 -Ofast

execute:
	./main.out


plotmc:
	python3 plotcycles.py
