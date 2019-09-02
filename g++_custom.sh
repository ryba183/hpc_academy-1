#!/usr/bin/env bash

input="${1}"
shift 1
g++ --pedantic -Wall -Wextra "${input}" -o "${input%.*}" "${@}"
