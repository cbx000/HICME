objects = build/testHIeB.o build/HIeB.o build/sq.o
CFLAGS = -std=c++0x -Wall -O
LDFLAGS = -lcuba -lgsl -lgslcblas -lm
BIN = bin/
BUILD = build/

bin/testHIeB: $(objects) $(BIN)
	g++ $(CFLAGS) -o $@ $(objects) $(LDFLAGS)
build/testHIeB.o: src/testHIeB.cpp src/HIeB.hpp $(BUILD) 
	g++ $(CFLAGS) -c src/testHIeB.cpp -o $@
build/HIeB.o: src/HIeB.cpp src/sq.h src/HIeB.hpp $(BUILD)
	g++ $(CFLAGS) -c src/HIeB.cpp -o $@
build/sq.o: src/sq.h src/sq.cpp $(BUILD)
	g++ $(CFLAGS) -c src/sq.cpp -o $@

$(BIN):
	mkdir -p $@
$(BUILD):
	mkdir -p $@
.PHONY: clean
clean:
	rm bin/testHIeB $(objects)
