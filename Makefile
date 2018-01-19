#objects = build/testHIeB.o build/HIeB.o build/sq.o
CFLAGS = -std=c++0x -Wall -O2
LDFLAGS = -lcuba -lgsl -lgslcblas -lm
BIN = bin/
BUILD = build/
all: bin/tabgen bin/testHIeB bin/eBgen
bin/eBgen: build/eBgen.o build/HIeB.o build/sq.o $(BIN)
	g++ $(CFLAGS) -o $@ build/eBgen.o build/HIeB.o build/sq.o $(LDFLAGS)
bin/tabgen: build/tabgen.o build/HIeB.o build/sq.o $(BIN)
	g++ $(CFLAGS) -o $@ build/tabgen.o build/HIeB.o build/sq.o $(LDFLAGS)
bin/testHIeB: build/testHIeB.o build/HIeB.o build/sq.o $(BIN)
	g++ $(CFLAGS) -o $@ build/testHIeB.o build/HIeB.o build/sq.o $(LDFLAGS)

build/eBgen.o: src/eBgen.cpp src/HIeB.hpp $(BUILD)
	g++ $(CFLAGS) -c src/eBgen.cpp -o $@
build/tabgen.o: src/tabgen.cpp src/HIeB.hpp $(BUILD) 
	g++ $(CFLAGS) -c src/tabgen.cpp -o $@
build/testHIeB.o: src/testHIeB.cpp src/HIeB.hpp $(BUILD) 
	g++ $(CFLAGS) -c src/testHIeB.cpp -o $@
build/HIeB.o: src/HIeB.cpp src/sq.hpp src/HIeB.hpp $(BUILD)
	g++ $(CFLAGS) -c src/HIeB.cpp -o $@
build/sq.o: src/sq.hpp src/sq.cpp $(BUILD)
	g++ $(CFLAGS) -c src/sq.cpp -o $@

$(BIN):
	mkdir -p $@
$(BUILD):
	mkdir -p $@
.PHONY: clean
clean:
	rm -f bin/*
	rm -f build/*.o
