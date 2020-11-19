all: compile execute

compile:
	g++ -o main.out $(wildcard *.cpp) -larmadillo -std=c++11 -Ofast -lomp -Xpreprocessor -fopenmp

execute:
	./main.out


plotmc:
	python3 plotcycles.py

hist:
	python3 plothist.py
