VPATH = src:build
objects = build/testHIeB.o build/HIeB.o build/sq.o
CFLAGS = -std=c++0x -Wall -O
LDFLAGS = -lcuba -lgsl -lgslcblas -lm

bin/testHIeB: $(objects)
	g++ $(CFLAGS) -o $@ $(objects) $(LDFLAGS)
build/testHIeB.o: testHIeB.cpp HIeB.hpp 
	g++ $(CFLAGS) -c src/testHIeB.cpp -o $@
build/HIeB.o: HIeB.cpp sq.h HIeB.hpp
	g++ $(CFLAGS) -c src/HIeB.cpp -o $@
build/sq.o: sq.h sq.cpp 
	g++ $(CFLAGS) -c src/sq.cpp -o $@

.PHONY: clean
clean:
	rm bin/testHIeB $(objects)
