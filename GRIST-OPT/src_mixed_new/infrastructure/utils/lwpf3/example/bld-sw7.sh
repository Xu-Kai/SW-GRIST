set -x
sw7gcc -mhost lwpf_example.c -c -o host.o -msimd -I.. &&
    sw7gcc -mslave lwpf_example.c -c -o slave.o -DCPE -msimd -O3 -g -I.. &&
    sw7gcc -mhybrid host.o slave.o