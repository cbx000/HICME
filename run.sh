#!/bin/bash
mkdir -p data
make
./bin/testHIeB
./bin/tabgen
