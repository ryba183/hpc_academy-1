#!/usr/bin/env bash
N=6
for i in $(seq 3 2 100) $(seq 101 10 1000) $(seq 1001 100 10000); do
 ((e=e%N)); ((e++==0)) && wait
 {
 mid=$(( (i / 2) + 1))
 vmid=$(./euler_equations_commandline -n "${i}" -t 0.15 -d 1 -D 1 -v -2 -V 2 -p 0.4 -P 0.4 | cut -f 2 | sed -n "$mid p")
 echo -e ${mid}"\t"${vmid}
 >&2 echo "${mid}"
} &
done; wait
#./euler_equations_commandline -n 1000 -t 0.15 -d 1 -D 1 -v -2 -V 2 -p 0.4 -P 0.4 | cut -f 2 | sed -n "500 p"

