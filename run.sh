#!/bin/bash
mkdir -p data
make
./bin/testHIeB
./bin/tabgen
./bin/eBgen
./bin/cme_ori
