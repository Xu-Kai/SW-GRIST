#include <omp.h>
#include <stdio.h>
int main(){
  double a = 0, b=0, c = 0;
  float d = 1;
  #pragma  omp parallel reduction(max:a) reduction(-:b) reduction(min:c) reduction(*:d)
  {
    a = c = d = omp_get_thread_num() + 1;
    d *= 0.1;
    b = -1;
  }
  printf("%f %f %f %f\n", a, b, c, d);
  #pragma  omp target parallel reduction(max:a) reduction(-:b) reduction(min:c) reduction(*:d)
  {
    a = c = d = omp_get_thread_num() + 1;
    d *= 0.1;
    b = -1;
  }
  printf("%f %f %f %f\n", a, b, c, d);
}
//{name = 0x4ba6ad "casw", opcode = 57344, mask = 4026593280, flags = 36864, operands = "\002\030\f\000"}
//{name = 0x4ba839 "addw", opcode = 1073741824, mask = 4026597120, flags = 4096, operands = "\002\003\004\000"}
//{name = 0x4ba6a3 "faal", opcode = 49152, mask = 4026593280, flags = 36864, operands = "\002\030\f\000"}