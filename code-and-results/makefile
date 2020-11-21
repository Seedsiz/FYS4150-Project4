all: compile execute

compile:
	g++-10 -Wall -o main.out $(wildcard *.cpp) -larmadillo -fopenmp -std=c++11 #-O2

execute:
	./main.out

plotmc:
	python3 plotcycles.py

hist:
	python3 plothist.py

expvalues:
	python3 plotexpvalues.py

time:
	python3 plottime.py

spin:
		python3 plotspin.py
