all: compile execute

compile:
	c++ -Wall -o main.out $(wildcard *.cpp) -larmadillo -lomp -std=c++11 -Ofast

execute:
	./main.out


plotmc:
	python3 plotcycles.py
