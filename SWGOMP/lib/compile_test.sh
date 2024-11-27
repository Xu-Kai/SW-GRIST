sw9gcc -mslave $1 -o $1.slave.o -c -DTEST -g -O3 &&
sw9gcc -mhost $1 -o $1.host.o -c -DTEST -g -O3 &&
sw9gcc -mslave tasktree.c -o tasktree.o -c -g -O3 &&
sw9gcc -mhybrid $1.host.o $1.slave.o tasktree.o -o $1.exe