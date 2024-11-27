set -x
sw5cc -host lwpf_example.c -c -o host.o -msimd -I.. &&
sw5cc -slave lwpf_example.c -c -o slave.o -msimd -O3 -I.. &&
sw5cc -hybrid host.o slave.o

