#!/usr/bin/env bash

g++ --pedantic -Wall -Wextra "${1}" -o "${1%.*}"
