set -x
sw9gcc -mhost lwpf_example.c -c -o host.o -msimd -I.. &&
    sw9gcc -mslave lwpf_example.c -c -o slave.o -DCPE -msimd -O3 -g -I.. &&
    sw9gcc -mhybrid host.o slave.o
