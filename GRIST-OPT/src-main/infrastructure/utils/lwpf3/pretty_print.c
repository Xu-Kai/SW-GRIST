#include <string.h>
#include <ctype.h>
#ifndef LWPF_FW
#define LWPF_FW 17
#endif
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

void pretty_print(char *to, char *s, int width, int max_lines) {
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
	  char *dots = "......";
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
#include <stdio.h>

int main(){
  char buf[LWPF_FW * 3];
  pretty_print(buf, "sumupibatchloop1masterckernel", LWPF_FW, 1);
  putchar('|');
  fwrite(buf, 1, LWPF_FW, stdout);
  putchar('|');
  puts("");

  pretty_print(buf, "sum_up_ibatch_loop1_master_c_kernel", LWPF_FW, 1);
  putchar('|');
  fwrite(buf, 1, LWPF_FW, stdout);
  putchar('|');
  puts("");
#define PCR_OPT_DEF(ev, val) {				\
    char *s = #ev;					\
    pretty_print(buf, s + 4, LWPF_FW, 3);		\
    putchar('|');					\
    fwrite(buf, 1, LWPF_FW, stdout);			\
    putchar('|');					\
    puts("");						\
    putchar('|');					\
    fwrite(buf + LWPF_FW, 1, LWPF_FW, stdout);		\
    putchar('|');					\
    puts("");						\
    putchar('|');					\
    fwrite(buf + LWPF_FW * 2, 1, LWPF_FW, stdout);	\
    putchar('|');					\
    puts("");						\
    puts("");						\
  }
  PCR_OPT_DEF(PC0_CYCLE                           , 0x00);

  PCR_OPT_DEF(PC1_INST                            , 0x00);
  PCR_OPT_DEF(PC1_NULL1                           , 0x01);
  PCR_OPT_DEF(PC1_INST_LAN                        , 0x02);
  PCR_OPT_DEF(PC1_NULL3                           , 0x03);
  PCR_OPT_DEF(PC1_NULL4                           , 0x04);
  PCR_OPT_DEF(PC1_NULL5                           , 0x05);
  PCR_OPT_DEF(PC1_INST_CBR                        , 0x06);
  PCR_OPT_DEF(PC1_INST_MEMACC                     , 0x07);
  PCR_OPT_DEF(PC1_INST_LDMACCESS                  , 0x08);
  PCR_OPT_DEF(PC1_DCACHE_HIT_LDM                  , 0x09);
  PCR_OPT_DEF(PC1_DCACHE_ACCESS                   , 0x0a);
  PCR_OPT_DEF(PC1_RMACCESS                        , 0x0b);
  PCR_OPT_DEF(PC1_INST_PIPE0                      , 0x0c);
  PCR_OPT_DEF(PC1_NULL13                          , 0x0d);
  PCR_OPT_DEF(PC1_INST_SCALAR_FLOAT_OTHER         , 0x0e);
  PCR_OPT_DEF(PC1_NULL15                          , 0x0f);
  PCR_OPT_DEF(PC1_REQ_DMA_RMA                     , 0x10);
  PCR_OPT_DEF(PC1_REQ_DMA                         , 0x11);
  PCR_OPT_DEF(PC1_REQ_RMA                         , 0x12);
  PCR_OPT_DEF(PC1_REQ_IO                          , 0x13);
  PCR_OPT_DEF(PC1_INST_SCALAR_HALF_FLOAT_OTHER    , 0x14);
  PCR_OPT_DEF(PC2_L1IC_ACCESS                     , 0x00);
  PCR_OPT_DEF(PC2_L1IC_MISSTIME                   , 0x01);
  PCR_OPT_DEF(PC2_CYC_DOUBLELAUNHC                , 0x02);
  PCR_OPT_DEF(PC2_NULL3                           , 0x03);
  PCR_OPT_DEF(PC2_NULL4                           , 0x04);
  PCR_OPT_DEF(PC2_INST_JUMP_CON                   , 0x05);
  PCR_OPT_DEF(PC2_PRE_SEQUENCE                    , 0x06);
  PCR_OPT_DEF(PC2_LDM_PIPEDIRECT                  , 0x07);
  PCR_OPT_DEF(PC2_CYC_LDMCONFLICT                 , 0x08);
  PCR_OPT_DEF(PC2_DCACHE_EVICT_LDM                , 0x09);
  PCR_OPT_DEF(PC2_DCACHE_MISS                     , 0x0a);
  PCR_OPT_DEF(PC2_GLD                             , 0x0b);
  PCR_OPT_DEF(PC2_INST_PIPE0_INT                  , 0x0c);
  PCR_OPT_DEF(PC2_INST_SCALAR_INT                 , 0x0d);
  PCR_OPT_DEF(PC2_INST_SCALAR_FLOAT_ADDETC        , 0x0e);
  PCR_OPT_DEF(PC2_REQ_RSYN                        , 0x0f);
  PCR_OPT_DEF(PC2_REQ_DMA                         , 0x10);
  PCR_OPT_DEF(PC2_REQ_DMA_PUT                     , 0x11);
  PCR_OPT_DEF(PC2_REQ_RMA_PUT                     , 0x12);
  PCR_OPT_DEF(PC2_REQ_L1IC_RAM_READ               , 0x13);
  PCR_OPT_DEF(PC2_INST_SCALAR_HALF_FLOAT_ADDETC   , 0x14);
  PCR_OPT_DEF(PC3_L1IC_MISS                       , 0x00);
  PCR_OPT_DEF(PC3_L0IC_BUBBLE                     , 0x01);
  PCR_OPT_DEF(PC3_CYC_LAUNHC_NONE                 , 0x02);
  PCR_OPT_DEF(PC3_NULL3                           , 0x03);
  PCR_OPT_DEF(PC3_NULL4                           , 0x04);
  PCR_OPT_DEF(PC3_FAIL_PREJUMP_INST_CON           , 0x05);
  PCR_OPT_DEF(PC3_PREJUMP_FORWARD                 , 0x06);
  PCR_OPT_DEF(PC3_DCACHE_HIT_LDM                  , 0x07);
  PCR_OPT_DEF(PC3_LDM_PIPE                        , 0x08);
  PCR_OPT_DEF(PC3_WRITELDM_DMA_RMA                , 0x09);
  PCR_OPT_DEF(PC3_REQ_DCACHE_EVICT                , 0x0a);
  PCR_OPT_DEF(PC3_GST                             , 0x0b);
  PCR_OPT_DEF(PC3_INST_PIPE0_FLOAT                , 0x0c);
  PCR_OPT_DEF(PC3_INST_SCALAR_FLOAT               , 0x0d);
  PCR_OPT_DEF(PC3_INST_SCALAR_FLOAT_MULETC        , 0x0e);
  PCR_OPT_DEF(PC3_REQ_CSYN                        , 0x0f);
  PCR_OPT_DEF(PC3_REQ_RMA                         , 0x10);
  PCR_OPT_DEF(PC3_REQ_DMA_GET                     , 0x11);
  PCR_OPT_DEF(PC3_REQ_RMA_GET                     , 0x12);
  PCR_OPT_DEF(PC3_REQ_LDM_RAM_READ                , 0x13);
  PCR_OPT_DEF(PC3_INST_SCALAR_HALF_FLOAT_MULETC   , 0x14);
  PCR_OPT_DEF(PC4_DCACHE_ACCESS                   , 0x00);
  PCR_OPT_DEF(PC4_CYC_LAUNHC_NONE_BUFFER          , 0x01);
  PCR_OPT_DEF(PC4_INST_PIPE0LAUNHC                , 0x02);
  PCR_OPT_DEF(PC4_CYC_PIPE0LAUNHCNONE_INST        , 0x03);
  PCR_OPT_DEF(PC4_L1IC_MISSTIME                   , 0x04);
  PCR_OPT_DEF(PC4_BR                              , 0x05);
  PCR_OPT_DEF(PC4_PREJUMP_BACK                    , 0x06);
  PCR_OPT_DEF(PC4_REQ_RSCATTER_MEMACC             , 0x07);
  PCR_OPT_DEF(PC4_CYC_LSCATTER_MEMACC_UNLDM       , 0x08);
  PCR_OPT_DEF(PC4_READLDM_DMA_RMA                 , 0x09);
  PCR_OPT_DEF(PC4_REQ_DCACHE_EVICT_ALONE          , 0x0a);
  PCR_OPT_DEF(PC4_RLD                             , 0x0b);
  PCR_OPT_DEF(PC4_INST_PIPE0_MEMACC               , 0x0c);
  PCR_OPT_DEF(PC4_INST_SCALAR_MEMACC              , 0x0d);
  PCR_OPT_DEF(PC4_INST_SCALAR_FLOAT_DIVETC        , 0x0e);
  PCR_OPT_DEF(PC4_REQ_SYNC_P2P                    , 0x0f);
  PCR_OPT_DEF(PC4_REQ_DMA_FENCE                   , 0x10);
  PCR_OPT_DEF(PC4_NULL17                          , 0x11);
  PCR_OPT_DEF(PC4_REQ_RMA_FENCE                   , 0x12);
  PCR_OPT_DEF(PC4_REQ_GPR_READ                    , 0x13);
  PCR_OPT_DEF(PC5_DCACHE_MISS                     , 0x00);
  PCR_OPT_DEF(PC5_CYC_SYN_WAIT                    , 0x01);
  PCR_OPT_DEF(PC5_INST_PIPE1LAUNHC                , 0x02);
  PCR_OPT_DEF(PC5_CYC_PIPE1LAUNHCNONE_INST        , 0x03);
  PCR_OPT_DEF(PC5_INST_L0IC_READ                  , 0x04);
  PCR_OPT_DEF(PC5_JMP                             , 0x05);
  PCR_OPT_DEF(PC5_FAIL_PRE_SEQUENCE               , 0x06);
  PCR_OPT_DEF(PC5_CYC_LSCATTER_MEMACC_UNARBITRATE , 0x07);
  PCR_OPT_DEF(PC5_DCACHE_ACCESS_LDM               , 0x08);
  PCR_OPT_DEF(PC5_LDM_AUTOADDONE                  , 0x09);
  PCR_OPT_DEF(PC5_REQ_DCACHE_PREFETCH             , 0x0a);
  PCR_OPT_DEF(PC5_RST                             , 0x0b);
  PCR_OPT_DEF(PC5_INST_PIPE1                      , 0x0c);
  PCR_OPT_DEF(PC5_INST_VECTOR_INT                 , 0x0d);
  PCR_OPT_DEF(PC5_INST_VECTOR_FLOAT_OTHER         , 0x0e);
  PCR_OPT_DEF(PC5_NULL15                          , 0x0f);
  PCR_OPT_DEF(PC5_REQ_RMA_FENCE                   , 0x10);
  PCR_OPT_DEF(PC5_REQ_DMA_FENCE_TOTAL             , 0x11);
  PCR_OPT_DEF(PC5_REQ_RMA_FENCE_TOTAL             , 0x12);
  PCR_OPT_DEF(PC5_REQ_L1IC_RAM_READ               , 0x13);
  PCR_OPT_DEF(PC5_INST_VECTOR_HALF_FLOAT_OTHER    , 0x14);
  PCR_OPT_DEF(PC6_CYC_LAUNHCNONE_BUFFER           , 0x00);
  PCR_OPT_DEF(PC6_CYC_MEMB_WAIT                   , 0x01);
  PCR_OPT_DEF(PC6_NULL2                           , 0x02);
  PCR_OPT_DEF(PC6_CYC_PIPE0_NOINST                , 0x03);
  PCR_OPT_DEF(PC6_L0IC_BUBBLE                     , 0x04);
  PCR_OPT_DEF(PC6_INST_CONFLAG                    , 0x05);
  PCR_OPT_DEF(PC6_FAIL_PREJUMP_FORWARD            , 0x06);
  PCR_OPT_DEF(PC6_CYC_DCACHE_CRQ                  , 0x07);
  PCR_OPT_DEF(PC6_LDM_TBOX                        , 0x08);
  PCR_OPT_DEF(PC6_ACCLDM_REMOTERLD                , 0x09);
  PCR_OPT_DEF(PC6_AVAILREQ_DCACHE_PREFETCH        , 0x0a);
  PCR_OPT_DEF(PC6_REQ_DCACHE_FILL                 , 0x0b);
  PCR_OPT_DEF(PC6_INST_PIPE1_INT                  , 0x0c);
  PCR_OPT_DEF(PC6_INST_VECTOR_FLOAT               , 0x0d);
  PCR_OPT_DEF(PC6_INST_VECTOR_FLOAT_ADDS          , 0x0e);
  PCR_OPT_DEF(PC6_REQ_OUTSYN                      , 0x0f);
  PCR_OPT_DEF(PC6_REQ_FENCE_TOTAL                 , 0x10);
  PCR_OPT_DEF(PC6_INST_BUFFIRSTDMARMA_0_0         , 0x11);
  PCR_OPT_DEF(PC6_NULL18                          , 0x12);
  PCR_OPT_DEF(PC6_REQ_LDM_RAM_READ                , 0x13);
  PCR_OPT_DEF(PC6_INST_VECTOR_HALF_FLOAT_ADDS     , 0x14);
  PCR_OPT_DEF(PC7_INST_L0IC_READ                  , 0x00);
  PCR_OPT_DEF(PC7_INST_BUFFIRSTDMARMA_0_0         , 0x01);
  PCR_OPT_DEF(PC7_CYC_INSTBUFFEREMPTY             , 0x02);
  PCR_OPT_DEF(PC7_CYC_PIPE1_NOINST                , 0x03);
  PCR_OPT_DEF(PC7_NULL4                           , 0x04);
  PCR_OPT_DEF(PC7_NULL5                           , 0x05);
  PCR_OPT_DEF(PC7_FAIL_PREJUMP_BACK               , 0x06);
  PCR_OPT_DEF(PC7_CYC_RSCATTER_MEMACC_UNARBITRATE , 0x07);
  PCR_OPT_DEF(PC7_LDM_RMW_WR                      , 0x08);
  PCR_OPT_DEF(PC7_ACCLDM_REMOTERST                , 0x09);
  PCR_OPT_DEF(PC7_CYC_DCACHE_CRQ                  , 0x0a);
  PCR_OPT_DEF(PC7_REQ_DCACHE_EVICT                , 0x0b);
  PCR_OPT_DEF(PC7_INST_PIPE1_MEMACC               , 0x0c);
  PCR_OPT_DEF(PC7_INST_VECTOR_MEMACC              , 0x0d);
  PCR_OPT_DEF(PC7_INST_VECTOR_FLOAT_MULS          , 0x0e);
  PCR_OPT_DEF(PC7_REQ_WRCH                        , 0x0f);
  PCR_OPT_DEF(PC7_REQ_FENCE_TOTAL_WAIT            , 0x10);
  PCR_OPT_DEF(PC7_NULL17                          , 0x11);
  PCR_OPT_DEF(PC7_NULL18                          , 0x12);
  PCR_OPT_DEF(PC7_REQ_GPR_READ                    , 0x13);
  PCR_OPT_DEF(PC7_INST_VECTOR_HALF_FLOAT_MULS     , 0x14);

}
