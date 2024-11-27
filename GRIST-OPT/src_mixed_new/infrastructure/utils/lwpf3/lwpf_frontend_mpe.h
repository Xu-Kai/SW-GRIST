#include <athread.h>
#include <stdio.h>
#include <ctype.h>

#define U(x) extern const char *lwpf_kernel_names_ ## x[];
LWPF_UNITS
#undef U

#define U(x) extern long lwpf_global_counter_ ## x[];
LWPF_UNITS
#undef U

#define U(x) extern long lwpf_kernel_count_ ## x;
LWPF_UNITS
#undef U
#define U(x) EXTERNC void slave_lwpf_init_ ## x(void *);
LWPF_UNITS
#undef U

#define U(x) evt_conf_t lwpf_evt_config_ ## x;
LWPF_UNITS;
#undef U
#ifndef LWPF_SPAWN
#define LWPF_SPAWN(func, arg) __real_athread_spawn((void*)(func), arg, 1)
#endif
#ifndef LWPF_JOIN
#define LWPF_JOIN athread_join
#endif

static inline void lwpf_init(void *conf){
#define U(x) LWPF_SPAWN(slave_lwpf_init_ ## x, conf); LWPF_JOIN();
  LWPF_UNITS
#undef U
}

/*
//GMK1
#ifndef LWPF_NOCOLOR
#define LWPF_KE   "\x1B[31mE\x1B[0m"
#define LWPF_KP   "\x1B[35mP\x1B[0m"
#define LWPF_KT   "\x1B[33mT\x1B[0m"
#define LWPF_KG   "\x1B[32mG\x1B[0m"
#define LWPF_KM   "\x1B[36mM\x1B[0m"
#define LWPF_KK   "\x1B[34mK\x1B[0m"
#define LWPF_K1   "\x1B[0m."
#define LWPF_KNRM "\x1B[0m"
#else
#define LWPF_KE   "T"
#define LWPF_KP"P"
#define LWPF_KK   "E"
#define LWPF_KG   "G"
#define LWPF_KM   "M"
#define LWPF_KK   "K"
#define LWPF_K1   "."
#define LWPF_KNRM ""
#endif
#define LWPF_DE   1000000000000000000L
#define LWPF_OE    500000000000000000L
#define LWPF_DP   1000000000000000L
#define LWPF_OP    500000000000000L
#define LWPF_DT   1000000000000L
#define LWPF_OT    500000000000L
#define LWPF_DG   1000000000L
#define LWPF_OG    500000000L
#define LWPF_DM   1000000L
#define LWPF_OM    500000L
#define LWPF_DK   1000L
#define LWPF_OK    500L
#define LWPF_D1   1L
#define LWPF_O1   0L
static char *LWPF_KN[] = {LWPF_K1, LWPF_KK, LWPF_KM, LWPF_KG, LWPF_KT, LWPF_KP, LWPF_KE};
static long LWPF_DN[] = {LWPF_D1, LWPF_DK, LWPF_DM, LWPF_DG, LWPF_DT, LWPF_DP, LWPF_DE};
static long LWPF_ON[] = {LWPF_O1, LWPF_OK, LWPF_OM, LWPF_OG, LWPF_OT, LWPF_OP, LWPF_OE};
static long LWPF_WN[] = {0, 3, 6, 9, 999};

static inline void print_unit_num_wx(FILE *output, long value, int w){
  int ptrw = 0;
  while (LWPF_WN[ptrw + 1] < w) ptrw ++;
  int ptrd = 0;
  while (ptrd < 7 && (value + LWPF_ON[ptrd]) / LWPF_DN[ptrd] >= LWPF_DN[ptrw] * 10) ptrd ++;
  fprintf(output, "%*lld%s%s", w - 1, (value + LWPF_ON[ptrd]) / LWPF_DN[ptrd], LWPF_KN[ptrd], LWPF_KNRM);
}
*/
#ifndef LWPF_NOCOLOR
#define LWPF_KE   "\x1B[31mG\x1B[0m"
#define LWPF_KP   "\x1B[35mG\x1B[0m"
#define LWPF_KT   "\x1B[33mG\x1B[0m"
#define LWPF_KG   "\x1B[32mG\x1B[0m"
#define LWPF_KM   "\x1B[36mM\x1B[0m"
#define LWPF_KK   "\x1B[34mK\x1B[0m"
#define LWPF_K1   "\x1B[0m."
#define LWPF_KNRM "\x1B[0m"
#else
#define LWPF_KE   "E"
#define LWPF_KP   "P"
#define LWPF_KT   "T"
#define LWPF_KG   "G"
#define LWPF_KM   "M"
#define LWPF_KK   "K"
#define LWPF_K1   "."
#define LWPF_KNRM ""
#endif
#define LWPF_DE   1000000000000000000
#define LWPF_OE    500000000000000000
#define LWPF_DP   1000000000000000
#define LWPF_OP    500000000000000
#define LWPF_DT   1000000000000
#define LWPF_OT    500000000000
#define LWPF_DG   1000000000
#define LWPF_OG    500000000
#define LWPF_DM   1000000
#define LWPF_OM    500000
#define LWPF_DK   1000
#define LWPF_OK    500
#define LWPF_D1   1
#define LWPF_O1   0
static const char *LWPF_KN[] = {LWPF_K1, LWPF_KK, LWPF_KM, LWPF_KG, LWPF_KT, LWPF_KP, LWPF_KE};
static long LWPF_DN[] =  {LWPF_D1, LWPF_DK, LWPF_DM, LWPF_DG, LWPF_DT, LWPF_DP, LWPF_DE};
static long LWPF_ON[] =  {LWPF_O1, LWPF_OK, LWPF_OM, LWPF_OG, LWPF_OT, LWPF_OP, LWPF_OE};
static long LWPF_WN[] = {0, 3, 6, 9, 999};

static inline void print_unit_num_wx(FILE *output, long value, int w){
  int ptrw = 0;
  while (LWPF_WN[ptrw + 1] < w) ptrw ++;
  int ptrd = 0;
  while (ptrd < 7 && (value + LWPF_ON[ptrd]) / LWPF_DN[ptrd] >= LWPF_DN[ptrw] * 10) ptrd ++;
  fprintf(output, "%*lld%s%s", w - 1, (value + LWPF_ON[ptrd]) / LWPF_DN[ptrd], LWPF_KN[ptrd], LWPF_KNRM);
}

void center_print(FILE *output, const char *s, int width) {
  int length = strlen(s);
  int i;
  for (i=0; i< (width-length) / 2; i++) {
    fputs(" ", output);
  }
  fputs(s, output);
  i += length;
  for (; i < width; i++) {
    fputs(" ", output);
  }
}

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

void make_center(char *to, int p, int width){
  int shift = (width - p) / 2;
  int k;
  for (k = p - 1; k >= 0; k --){
    to[k + shift] = to[k];
  }
  for (k = 0; k < shift; k ++)
    to[k] = ' ';
}

void make_vcenter(char *to, int l, int width, int max_lines){
  int shift = (max_lines - l) / 2;
  if (shift == 0) return;
  int k;
  for (k = l - 1; k >= 0; k --){
    memcpy(to + (k + shift) * width, to + k * width, width * sizeof(char));
  }
  for (k = 0; k < shift; k ++)
    memset(to + k * width, ' ', width * sizeof(char));
}

/*void pretty_print(char *to, char *s, int width, int max_lines) {
  memset(to, ' ', width * max_lines * sizeof(char));
  int i = 0;
  int p = 0;
  int l = 0;
  while (1){
    while (s[i] == '_') i ++;
    if (!s[i] || s[i] == '\n') break;
    int j = 0;
    while (isalnum(s[i + j])) j ++;

    int needed = p == 0 ? j : j + 1;
    if (p + needed < width) {
      p ++;
      strncpy(to + l * width + p, s + i, j);
      p += j;
      i += j;
    } else {
      make_center(to + l * width, p, width);
      l ++;
      p = 0;
    }
  }
  if (p == 0)
    l --;
  else
    make_center(to + l * width, p, width);
  l ++;
  make_vcenter(to, l, width, max_lines);
}
*/
void pretty_print(char *to, const char *s, int width, int max_lines) {
  memset(to, ' ', width * max_lines * sizeof(char));
  int i = 0;
  int p = 0;
  int l = 0;
  while (1){
    while (s[i] == '_') i ++;
    if (!s[i] || s[i] == '\n') break;
    int j = 0;
    while (isalnum(s[i + j])) j ++;

    int needed = p == 0 ? j : j + 1;
    if (p + needed < width) {
      if (p > 0)
	p ++;
      memcpy(to + l * width + p, s + i, j);
      p += j;
      i += j;
    } else {
      if (p > 0){
	if (l == max_lines - 1) {
	  const char *dots = "......";
	  if (p < width){
	    if (width - p < 6){
	      memcpy(to + l * width + p, dots, width - p);
	      p = width;
	    } else {
	      memcpy(to + l * width + p, dots, 6);
	      p += 6;
	    }
	  }
	}
	make_center(to + l * width, p, width);
	l ++;
	p = 0;
	if (l >= max_lines) break;
      }
      if (j > width) {
	int nprint = (width - 3) / 2;
	memcpy(to + l * width + p, s + i, nprint);
	memcpy(to + l * width + p + nprint, "...", 3);
	memcpy(to + l * width + p + nprint + 3, s + i + j - nprint, nprint);
	p = width;
	i += width;
      }
    }
  }
  if (p == 0)
    l --;
  else
    make_center(to + l * width, p, width);
  l ++;
  make_vcenter(to, l, width, max_lines);
}

#ifndef LWPF_FW
#define LWPF_FW 17
#endif
static inline long get_max(long *in, int stride, int count){
  int i;
  long ret = 0;
  for (i = 0; i < count; i ++){
    ret = max(ret, in[i * stride]);
  }
  return ret;
}
static inline long get_avg(long *in, int stride, int count){
  int i;
  long ret = 0;
  for (i = 0; i < count; i ++){
    ret = ret + in[i * stride];
  }
  return (ret + count / 2) / count;
}
static inline long get_min(long *in, int stride, int count){
  int i;
  long ret = 0x7fffffffffffffffL;
  for (i = 0; i < count; i ++){
    ret = min(ret, in[i * stride]);
  }
  return ret;
}
static inline void print_hline(FILE *output, int nfield, int width){
  int j, k;
  fputc('+', output);
  for (j = 0; j < nfield; j ++){
    for (k = 0; k < width; k ++)
      fputc('-', output);
    fputc('+', output);
  }
  fputs("\n", output);
}
/*
static inline void lwpf_report_summary_one(FILE *output, char *uname, long *lwpf_global_counter, char **lwpf_kernel_names, long lwpf_kernel_count, evt_conf_t *conf){
  printf("one: %s %p\n", uname, lwpf_global_counter);
  fprintf(output, "LWPF kernel summary for unit %s\n", uname);
  int i, j, k, l;
  int npcr = 0;
  int pcrs[NPCR];
  char ***pcr_names = PCR_NAMES;
  for (j = 0; j < NPCR; j ++)
    if (conf->pc_mask & (1 << j))
      pcrs[npcr ++] = j;
  char hdr_out[1 + NPCR][LWPF_FW * 3];
  pretty_print(hdr_out[0], "KERNEL", LWPF_FW, 3);
  for (j = 0; j < npcr; j ++){
    pretty_print(hdr_out[j + 1], pcr_names[pcrs[j]][conf->evt[pcrs[j]]] + 4, LWPF_FW, 3);
  }

  print_hline(output, npcr + 1, LWPF_FW);
  
  for (l = 0; l < 3; l ++){
    fprintf(output, "|");
    for (j = 0; j < npcr + 1; j ++){
      fwrite(hdr_out[j] + LWPF_FW * l, 1, LWPF_FW, output);
      fprintf(output, "|");
    }
    fprintf(output, "\n");
  }
  print_hline(output, npcr + 1, LWPF_FW);
  char ker_out[LWPF_FW * 10];
  int num_width = (LWPF_FW - 2) / 3;
  for (i = 0; i < lwpf_kernel_count; i ++){
    //pretty_print(ker_out, lwpf_kernel_names[i], LWPF_FW, 1);

    puts(lwpf_kernel_names[i]);

    //fprintf(output, "|");
   // puts("wyj_lwpf kernel print1");
    //fwrite(ker_out, 1, LWPF_FW, output);
   // puts("wyj_lwpf kernel print2");
    //fprintf(output, "|");
   // puts("wyj_lwpf kernel print3");
    printf("%p\n", lwpf_global_counter);
    for (j = 0; j < npcr; j ++){
   // puts("wyj_lwpf kernel print4");
      long max = get_max(&lwpf_global_counter[i * NPCR + pcrs[j]], lwpf_kernel_count * NPCR, 64);
   // puts("wyj_lwpf kernel print5");
      long min = get_min(&lwpf_global_counter[i * NPCR + pcrs[j]], lwpf_kernel_count * NPCR, 64);
      long avg = get_avg(&lwpf_global_counter[i * NPCR + pcrs[j]], lwpf_kernel_count * NPCR, 64);
      #define e2 2000000000.0
      printf("%-22s:", pcr_names[pcrs[j]][conf->evt[pcrs[j]]] + 4);
      printf("max=%20lld, min=%20lld, avg=%20lld\n", max, min, avg);
      //printf("max=%lf, min=%lf, avg=%lf\n", max/e2, min/e2, avg/e2);
     // print_unit_num_wx(output, avg, num_width);
     // fprintf(output, "|");
     // print_unit_num_wx(output, min, num_width);
     // fprintf(output, "|");
     // print_unit_num_wx(output, max, num_width);
     // fprintf(output, "|");
    }
    //puts("wyj_lwpf kernel print before fputs");
    //fputs("\n", output);
    //puts("wyj_lwpf kernel print end");
  }
}
*/
static inline void lwpf_report_summary_one(FILE *output, const char *uname, long *lwpf_global_counter, const char **lwpf_kernel_names, long lwpf_kernel_count, evt_conf_t *conf){
  fprintf(output, "LWPF kernel summary for unit %s\n", uname);
  int i, j, k, l;
  int npcr = 0;
  int pcrs[NPCR];
  const char ***pcr_names = PCR_NAMES;
  for (j = 0; j < NPCR; j ++)
    if (conf->pc_mask & (1 << j))
      pcrs[npcr ++] = j;

  char hdr_out[1 + NPCR][LWPF_FW * 3];
  pretty_print(hdr_out[0], "KERNEL", LWPF_FW, 3);
  for (j = 0; j < npcr; j ++){
    pretty_print(hdr_out[j + 1], pcr_names[pcrs[j]][conf->evt[pcrs[j]]] + 4, LWPF_FW, 3);
  }

  print_hline(output, npcr + 1, LWPF_FW);
  for (l = 0; l < 3; l ++){
    fprintf(output, "|");
    for (j = 0; j < npcr + 1; j ++){
      fwrite(hdr_out[j] + LWPF_FW * l, 1, LWPF_FW, output);
      fprintf(output, "|");
    }
    fprintf(output, "\n");
  }
  print_hline(output, npcr + 1, LWPF_FW);
  char ker_out[LWPF_FW];
  int num_width = (LWPF_FW - 2) / 3;
  for (i = 0; i < lwpf_kernel_count; i ++){
    pretty_print(ker_out, lwpf_kernel_names[i], LWPF_FW, 1);
    fprintf(output, "|");
    fwrite(ker_out, 1, LWPF_FW, output);
    fprintf(output, "|");
    for (j = 0; j < npcr; j ++){
      long max = get_max(&lwpf_global_counter[i * NPCR + pcrs[j]], lwpf_kernel_count * NPCR, 64);
      long min = get_min(&lwpf_global_counter[i * NPCR + pcrs[j]], lwpf_kernel_count * NPCR, 64);
      long avg = get_avg(&lwpf_global_counter[i * NPCR + pcrs[j]], lwpf_kernel_count * NPCR, 64);
      print_unit_num_wx(output, avg, num_width);
      fprintf(output, "|");
      print_unit_num_wx(output, min, num_width);
      fprintf(output, "|");
      print_unit_num_wx(output, max, num_width);
      fprintf(output, "|");
    }
    fputs("\n", output);
  }
}
static inline void lwpf_report_summary(FILE *output){
#define U(x) printf("summary: %s %p\n", #x, lwpf_global_counter_ ## x);
  LWPF_UNITS
#undef U

#define U(x) lwpf_report_summary_one(output, #x, lwpf_global_counter_ ## x, lwpf_kernel_names_ ## x, lwpf_kernel_count_ ## x, &lwpf_evt_config_ ## x);
  LWPF_UNITS
#undef U
}
					    
static inline void lwpf_report_detail_one(FILE *output, const char *uname, long *lwpf_global_counter, const char **lwpf_kernel_names, long lwpf_kernel_count, evt_conf_t *conf){
  fprintf(output, "LWPF kernel summary for unit %s\n", uname);
  int i, j, k, l;
  int npcr = 0;
  int pcrs[NPCR];
  const char ***pcr_names = PCR_NAMES;
  for (j = 0; j < NPCR; j ++)
    if (conf->pc_mask & (1 << j))
      pcrs[npcr ++] = j;

  char hdr_out[1 + NPCR][LWPF_FW * 3];
  pretty_print(hdr_out[0], "KERNEL", LWPF_FW, 3);
  for (j = 0; j < npcr; j ++){
    pretty_print(hdr_out[j + 1], pcr_names[pcrs[j]][conf->evt[pcrs[j]]] + 4, LWPF_FW, 3);
  }

  print_hline(output, npcr + 1, LWPF_FW);
  for (l = 0; l < 3; l ++){
    fprintf(output, "|");
    for (j = 0; j < npcr + 1; j ++){
      fwrite(hdr_out[j] + LWPF_FW * l, 1, LWPF_FW, output);
      fprintf(output, "|");
    }
    fprintf(output, "\n");
  }
  print_hline(output, npcr + 1, LWPF_FW);
  char ker_out[LWPF_FW * 3];
  int num_width = LWPF_FW;
  for (i = 0; i < lwpf_kernel_count; i ++){
    for (k = 0; k < 64; k ++){
      char ker_tmp[LWPF_FW * 3];
      sprintf(ker_tmp, "%s%d",lwpf_kernel_names[i], k);
      pretty_print(ker_out, ker_tmp, LWPF_FW, 1);
      fprintf(output, "|");
      fwrite(ker_out, 1, LWPF_FW, output);
      fprintf(output, "|");
      for (j = 0; j < npcr; j ++){
      	print_unit_num_wx(output, lwpf_global_counter[(k * lwpf_kernel_count + i) * NPCR + pcrs[j]], LWPF_FW);
      	fprintf(output, "|");
      }
      fputs("\n", output);
    }
    print_hline(output, npcr + 1, LWPF_FW);
  }
}

static inline void lwpf_report_detail(FILE *output){
#define U(x) lwpf_report_detail_one(output, #x, lwpf_global_counter_ ## x, lwpf_kernel_names_ ## x, lwpf_kernel_count_ ## x, &lwpf_evt_config_ ## x);
  LWPF_UNITS
#undef U
}
#undef max
#undef min
