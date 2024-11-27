#include <stdio.h>

int be2le(int in){
  union ib {
    int i;
    char b[4];
  } i, o;
  i.i = in;
  o.b[0] = i.b[3];
  o.b[1] = i.b[2];
  o.b[2] = i.b[1];
  o.b[3] = i.b[0];
  return o.i;
}
void get_bin_range_(void *buf, char *fp, long *displs, long *count){
  return get_bin_intbe_(buf, fp, displs, count);
}
void get_bin_intbe_(int *buf, char *fp, long *displs, long *size){
  FILE *f = fopen(fp, "rb");
  int count = *size / 4;
  fseek(f, *displs, SEEK_SET);
  fread(buf, 4, count, f);
  for (int i = 0; i < count; i ++) {
    buf[i] = be2le(buf[i]);
  }
  fclose(f);
}
