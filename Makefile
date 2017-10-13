testHIeB: testHIeB.o HIeB.o HIeB.hpp HIeB.cpp sq.c sq.h
	g++ -std=c++0x -Wall -O -o testHIeB testHIeB.o HIeB.o -lcuba -lgsl -lgslcblas -lm
testHIeB.o: testHIeB.cpp HIeB.hpp 
	g++ -std=c++0x -Wall -O -c testHIeB.cpp
HIeB.o: HIeB.cpp sq.h HIeB.hpp
	g++ -std=c++0x -Wall -O -c HIeB.cpp
sq.o: sq.h sq.c 
	g++ -std=c++0x -Wall -O -c sq.c

clean:
	rm *.o
	rm testHIeB
