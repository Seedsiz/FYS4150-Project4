all: compile execute

compile:
	clang++ -o main.out $(wildcard *.cpp) -larmadillo -std=c++11 -Ofast -lomp -Xpreprocessor

execute:
	./main.out


plotmc:
	python3 plotcycles.py

hist:
	python3 plothist.py
