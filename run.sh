#!/bin/bash
mkdir -p data/eta
make
./bin/testHIeB
./bin/tabgen
./bin/eBgen
./bin/cme_ori
./bin/cme_ratio
./bin/cme_eta
