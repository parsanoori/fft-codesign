#!/bin/bash

set -e

clear

echo "--> compiling sw"
arm-linux-gcc -static -Wall -std=c99 fft.c -o fft -O3 -lm

echo "--> run gplatform"
gplatform fft.fdl 2>/dev/null

echo "--> done"

