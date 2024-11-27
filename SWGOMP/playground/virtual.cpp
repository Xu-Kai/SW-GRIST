#include <stdio.h>
struct A {
#pragma omp declare target
  virtual void exec();
};
#pragma omp end declare target
void A::exec(){
  puts("A");
}
struct B : A{
  virtual void exec() override {
    puts("B");
  }
};
#pragma omp declare target
// #pragma omp declare target to(exec)
void exec(A &a){
  a.exec();
}
#pragma omp end declare target
int main(){
  A a;
  B b;
  auto G=[](){puts("G");};
  struct C : A{
  virtual void exec() override {
      puts("C");
    }
  };
  #pragma omp target firstprivate(a, b)
  {
    exec(a);
    exec(b);
  }
}