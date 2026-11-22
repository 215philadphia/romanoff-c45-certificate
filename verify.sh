#!/bin/sh
set -e
g++ -O2 -o seed seed.cpp
g++ -O2 -o mult mult.cpp
./seed 3 5 7 11 13 17 19 29 31 37 41 61 73 257 > seedY.txt   # ~3 min
python3 prove2.py clusters_113.json 100000000                # ~3 min, ~3 GB RAM
