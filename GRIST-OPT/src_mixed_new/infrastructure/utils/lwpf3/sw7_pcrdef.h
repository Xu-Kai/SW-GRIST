PCR_DEF_BEGIN(PC0)
//周期数:
PCR_OPT_DEF(PC0_CYCLE                          , 0x0)
PCR_DEF_END(PC0)

PCR_DEF_BEGIN(PC1)
//执行指令数:
PCR_OPT_DEF(PC1_INST                           , 0x0)
PCR_OPT_DEF(PC1_NULL1                          , 0x1)
//发射指令数:
PCR_OPT_DEF(PC1_INST_LAN                       , 0x2)
PCR_OPT_DEF(PC1_NULL3                          , 0x3)
PCR_OPT_DEF(PC1_NULL4                          , 0x4)
PCR_OPT_DEF(PC1_NULL5                          , 0x5)
//条件转移指令执行次数:
PCR_OPT_DEF(PC1_INST_CBR                       , 0x6)
//访存指令数(不含VLENMAS):
PCR_OPT_DEF(PC1_INST_MEMACC                    , 0x7)
//LDM访问总次数:
PCR_OPT_DEF(PC1_INST_LDMACCESS                 , 0x8)
PCR_OPT_DEF(PC1_NULL9                          , 0x9)
PCR_OPT_DEF(PC1_NULL10                         , 0xa)
//远程访问总次数:
PCR_OPT_DEF(PC1_RMACCESS                       , 0xb)
//流水线0执行指令数:
PCR_OPT_DEF(PC1_INST_PIPE0                     , 0xc)
PCR_OPT_DEF(PC1_NULL13                         , 0xd)
PCR_OPT_DEF(PC1_NULL14                         , 0xe)
PCR_OPT_DEF(PC1_NULL15                         , 0xf)
//DMA、RMA请求总数:
PCR_OPT_DEF(PC1_REQ_DMA_RMA                    , 0x10)
//DMA_PUT、DMA_GET请求数:
PCR_OPT_DEF(PC1_REQ_DMA                        , 0x11)
//RMA_PUT、RMA_GET 请求数:
PCR_OPT_DEF(PC1_REQ_RMA                        , 0x12)
//从核收到的IO请求总数(含对RAM的访问):
PCR_OPT_DEF(PC1_REQ_IO                         , 0x13)
PCR_DEF_END(PC1)

PCR_DEF_BEGIN(PC2)
//L1IC访问次数:
PCR_OPT_DEF(PC2_L1IC_ACCESS                    , 0x0)
//L1IC 脱靶总时间:
PCR_OPT_DEF(PC2_L1IC_MISSTIME                  , 0x1)
//双发射拍数:
PCR_OPT_DEF(PC2_CYC_DOUBLELAUNHC               , 0x2)
PCR_OPT_DEF(PC2_NULL3                          , 0x3)
PCR_OPT_DEF(PC2_NULL4                          , 0x4)
//条件转移指令执行次数:
PCR_OPT_DEF(PC2_INST_JUMP_CON                  , 0x5)
//顺序执行预测次数:
PCR_OPT_DEF(PC2_PRE_SEQUENCE                   , 0x6)
//流水线直接访问LDM的次数(含VLENMAS,不含RMW的写操作):
PCR_OPT_DEF(PC2_LDM_PIPEDIRECT                 , 0x7)
//LDM访问冲突总拍数:
PCR_OPT_DEF(PC2_CYC_LDMCONFLICT                , 0x8)
PCR_OPT_DEF(PC2_NULL9                          , 0x9)
PCR_OPT_DEF(PC2_NULL10                         , 0xa)
//GLD访问主存的次数:
PCR_OPT_DEF(PC2_GLD                            , 0xb)
//流水线0执行整数运算指令数:
PCR_OPT_DEF(PC2_INST_PIPE0_INT                 , 0xc)
//执行标量整数运算指令数:
PCR_OPT_DEF(PC2_INST_SCALAR_INT                , 0xd)
//标量浮点加减乘指令数:
PCR_OPT_DEF(PC2_INST_SCALAR_FLOAT_ADDETC       , 0xe)
//行同步请求数:
PCR_OPT_DEF(PC2_REQ_RSYN                       , 0xf)
//DMA_PUT、DMA_GET请求数:
PCR_OPT_DEF(PC2_REQ_DMA                        , 0x10)
//DMA_PUT请求数:
PCR_OPT_DEF(PC2_REQ_DMA_PUT                    , 0x11)
//RMA_PUT请求数:
PCR_OPT_DEF(PC2_REQ_RMA_PUT                    , 0x12)
//从核收到的对 L1IC_RAM 的 IO 写请求总数:
PCR_OPT_DEF(PC2_REQ_L1IC_RAM_READ              , 0x13)
PCR_DEF_END(PC2)

PCR_DEF_BEGIN(PC3)
//L1IC 脱靶次数:
PCR_OPT_DEF(PC3_L1IC_MISS                      , 0x0)
//L0IC 性能气泡:
PCR_OPT_DEF(PC3_L0IC_BUBBLE                    , 0x1)
//零发射拍数:
PCR_OPT_DEF(PC3_CYC_LAUNHC_NONE                , 0x2)
PCR_OPT_DEF(PC3_NULL3                          , 0x3)
PCR_OPT_DEF(PC3_NULL4                          , 0x4)
//条件转移预测失败次数:
PCR_OPT_DEF(PC3_FAIL_PREJUMP_INST_CON          , 0x5)
//前跳预测次数:
PCR_OPT_DEF(PC3_PREJUMP_FORWARD                , 0x6)
PCR_OPT_DEF(PC3_NULL7                          , 0x7)
//流水线访问 LDM 次数:
PCR_OPT_DEF(PC3_LDM_PIPE                       , 0x8)
//DMA/RMA 写访问 LDM 的次数:
PCR_OPT_DEF(PC3_WRITELDM_DMA_RMA               , 0x9)
PCR_OPT_DEF(PC3_NULL10                         , 0xa)
//GST 访问主存的次数:
PCR_OPT_DEF(PC3_GST                            , 0xb)
//流水线 0 执行浮点运算指令数:
PCR_OPT_DEF(PC3_INST_PIPE0_FLOAT               , 0xc)
//执行标量浮点运算指令数:
PCR_OPT_DEF(PC3_INST_SCALAR_FLOAT              , 0xd)
//标量浮点乘加指令数:
PCR_OPT_DEF(PC3_INST_SCALAR_FLOAT_MULETC       , 0xe)
//列同步请求数:
PCR_OPT_DEF(PC3_REQ_CSYN                       , 0xf)
//RMA_PUT、RMA_GET 请求数:
PCR_OPT_DEF(PC3_REQ_RMA                        , 0x10)
//DMA_GET 请求数:
PCR_OPT_DEF(PC3_REQ_DMA_GET                    , 0x11)
//RMA_GET 请求数:
PCR_OPT_DEF(PC3_REQ_RMA_GET                    , 0x12)
//从核收到的对 LDM_RAM 的 IO 写请求总数:
PCR_OPT_DEF(PC3_REQ_LDM_RAM_READ               , 0x13)
PCR_DEF_END(PC3)

PCR_DEF_BEGIN(PC4)
PCR_OPT_DEF(PC4_NULL0                          , 0x0)
//指令发射缓冲有指令但零发射拍数:
PCR_OPT_DEF(PC4_CYC_LAUNHC_NONE_BUFFER         , 0x1)
//流水线 0 发射指令数:
PCR_OPT_DEF(PC4_INST_PIPE0LAUNHC               , 0x2)
//流水线 0 有指令,但无法发射拍数:
PCR_OPT_DEF(PC4_CYC_PIPE0LAUNHCNONE_INST       , 0x3)
//L1IC 脱靶总时间:
PCR_OPT_DEF(PC4_L1IC_MISSTIME                  , 0x4)
//BR 执行次数:
PCR_OPT_DEF(PC4_BR                             , 0x5)
//后跳预测次数:
PCR_OPT_DEF(PC4_PREJUMP_BACK                   , 0x6)
//远程离散访存(GLD/GST/原子操作/RLD/RST)请求总数:
PCR_OPT_DEF(PC4_REQ_RSCATTER_MEMACC            , 0x7)
//本地离散访存(LLD/LST,不含FPU_ST)请求未被 LDM 仲裁上的拍数:
PCR_OPT_DEF(PC4_CYC_LSCATTER_MEMACC_UNLDM      , 0x8)
//DMA/RMA 读访问 LDM 的次数:
PCR_OPT_DEF(PC4_READLDM_DMA_RMA                , 0x9)
PCR_OPT_DEF(PC4_NULL10                         , 0xa)
//RLD 访问远程 LDM 的次数:
PCR_OPT_DEF(PC4_RLD                            , 0xb)
//流水线 0 执行访存指令数(不含VLENMAS):
PCR_OPT_DEF(PC4_INST_PIPE0_MEMACC              , 0xc)
//执行标量访存指令数:
PCR_OPT_DEF(PC4_INST_SCALAR_MEMACC             , 0xd)
//标量浮点除法平方根指令数:
PCR_OPT_DEF(PC4_INST_SCALAR_FLOAT_DIVETC       , 0xe)
//点对点同步请求:
PCR_OPT_DEF(PC4_REQ_SYNC_P2P                   , 0xf)
//DMA_栏栅请求数:
PCR_OPT_DEF(PC4_REQ_DMA_FENCE                  , 0x10)
PCR_OPT_DEF(PC4_NULL17                         , 0x11)
//RMA_栏栅请求数:
PCR_OPT_DEF(PC4_REQ_RMA_FENCE                  , 0x12)
//从核收到的对 GPR 的 IO 写请求总数:
PCR_OPT_DEF(PC4_REQ_GPR_READ                   , 0x13)
PCR_DEF_END(PC4)

PCR_DEF_BEGIN(PC5)
PCR_OPT_DEF(PC5_NULL0                          , 0x0)
//同步等待总拍数:
PCR_OPT_DEF(PC5_CYC_SYN_WAIT                   , 0x1)
//流水线 1 发射指令数:
PCR_OPT_DEF(PC5_INST_PIPE1LAUNHC               , 0x2)
//L0IC 读取指令数:
PCR_OPT_DEF(PC5_CYC_PIPE1LAUNHCNONE_INST       , 0x3)
//L0IC 读取指令数:
PCR_OPT_DEF(PC5_INST_L0IC_READ                 , 0x4)
//顺序执行预测失败次数:
PCR_OPT_DEF(PC5_JMP                            , 0x5)
//流水线 1 有指令,但无法发射拍数:
PCR_OPT_DEF(PC5_FAIL_PRE_SEQUENCE              , 0x6)
//本地离散访存(LLD/LST,不含FPU_ST)请求未被 LDM 仲裁上的拍数:
PCR_OPT_DEF(PC5_CYC_LSCATTER_MEMACC_UNARBITRATE, 0x7)
PCR_OPT_DEF(PC5_NULL8                          , 0x8)
//回答字自增 1 访问 LDM 的次数:
PCR_OPT_DEF(PC5_LDM_AUTOADDONE                 , 0x9)
PCR_OPT_DEF(PC5_NULL10                         , 0xa)
//RST 访问远程 LDM 的次数:
PCR_OPT_DEF(PC5_RST                            , 0xb)
//流水线 1 执行指令数:
PCR_OPT_DEF(PC5_INST_PIPE1                     , 0xc)
//执行向量整数运算指令数:
PCR_OPT_DEF(PC5_INST_VECTOR_INT                , 0xd)
//标量浮点转换指令数:
PCR_OPT_DEF(PC5_INST_SCALAR_FLOAT_TRANS        , 0xe)
PCR_OPT_DEF(PC5_NULL15                         , 0xf)
//RMA_栏栅请求数:
PCR_OPT_DEF(PC5_REQ_RMA_FENCE                  , 0x10)
//DMA_全栏栅请求数:
PCR_OPT_DEF(PC5_REQ_DMA_FENCE_TOTAL            , 0x11)
//RMA_全栏栅请求数:
PCR_OPT_DEF(PC5_REQ_RMA_FENCE_TOTAL            , 0x12)
//从核收到的对 L1IC_RAM 的 IO 写请求总数:
PCR_OPT_DEF(PC5_REQ_L1IC_RAM_READ              , 0x13)
PCR_DEF_END(PC5)

PCR_DEF_BEGIN(PC6)
//指令发射缓冲有指令但零发射拍数:
PCR_OPT_DEF(PC6_CYC_LAUNHCNONE_BUFFER          , 0x0)
//MEMB 等待总拍数:
PCR_OPT_DEF(PC6_CYC_MEMB_WAIT                  , 0x1)
PCR_OPT_DEF(PC6_NULL2                          , 0x2)
//流水线 0 无指令拍数:
PCR_OPT_DEF(PC6_CYC_PIPE0_NOINST               , 0x3)
//L0IC 性能气泡:
PCR_OPT_DEF(PC6_L0IC_BUBBLE                    , 0x4)
//带条件标记的运算指令执行次数:
PCR_OPT_DEF(PC6_INST_CONFLAG                   , 0x5)
//前跳预测失败次数:
PCR_OPT_DEF(PC6_FAIL_PREJUMP_FORWARD           , 0x6)
PCR_OPT_DEF(PC6_NULL7                          , 0x7)
//TBOX 访问 LDM 次数(DMA,RMA,回答字,外部 RLD/RST):
PCR_OPT_DEF(PC6_LDM_TBOX                       , 0x8)
//外部 RLD 请求读访问本地 LDM 的次数:
PCR_OPT_DEF(PC6_ACCLDM_REMOTERLD               , 0x9)
PCR_OPT_DEF(PC6_NULL10                         , 0xa)
PCR_OPT_DEF(PC6_NULL11                         , 0xb)
//流水线 1 执行整数运算指令数:
PCR_OPT_DEF(PC6_INST_PIPE1_INT                 , 0xc)
//执行向量浮点运算指令数(含VLENMAS):
PCR_OPT_DEF(PC6_INST_VECTOR_FLOAT              , 0xd)
//向量浮点加减乘指令数:
PCR_OPT_DEF(PC6_INST_VECTOR_FLOAT_ADDS         , 0xe)
//外部同步请求个数:
PCR_OPT_DEF(PC6_REQ_OUTSYN                     , 0xf)
//全栏栅请求数:
PCR_OPT_DEF(PC6_REQ_FENCE_TOTAL                , 0x10)
//指令发射缓冲首指令为访存、DMA或 RMA 指令,但对应信用为 0,且当拍零发射拍数:
PCR_OPT_DEF(PC6_INST_BUFFIRSTDMARMA_0_0        , 0x11)
PCR_OPT_DEF(PC6_NULL18                         , 0x12)
//从核收到的对 LDM_RAM 的 IO 写请求总数:
PCR_OPT_DEF(PC6_REQ_LDM_RAM_READ               , 0x13)
PCR_DEF_END(PC6)

PCR_DEF_BEGIN(PC7)
//L0IC 读取指令数:
PCR_OPT_DEF(PC7_INST_L0IC_READ                 , 0x0)
//指令发射缓冲首指令为访存、DMA或 RMA 指令,但对应信用为 0,且当拍零发射拍数:
PCR_OPT_DEF(PC7_INST_BUFFIRSTDMARMA_0_0        , 0x1)
//RMW_WR 访问 LDM 的次数:
PCR_OPT_DEF(PC7_CYC_INSTBUFFEREMPTY            , 0x2)
//指令发射缓冲无指令拍数:
PCR_OPT_DEF(PC7_CYC_PIPE1_NOINST               , 0x3)
PCR_OPT_DEF(PC7_NULL4                          , 0x4)
PCR_OPT_DEF(PC7_NULL5                          , 0x5)
//流水线 1 无指令拍数:
PCR_OPT_DEF(PC7_FAIL_PREJUMP_BACK              , 0x6)
//后跳预测失败次数:
PCR_OPT_DEF(PC7_CYC_RSCATTER_MEMACC_UNARBITRATE, 0x7)
//远程离散访存(GLD/GST/原子操作/RLD/RST)请求被 GGQ 满或未仲裁上阻塞的拍数:
PCR_OPT_DEF(PC7_LDM_RMW_WR                     , 0x8)
//外部 RST 请求写访问本地 LDM 的次数:
PCR_OPT_DEF(PC7_ACCLDM_REMOTERST               , 0x9)
PCR_OPT_DEF(PC7_NULL10                         , 0xa)
PCR_OPT_DEF(PC7_NULL11                         , 0xb)
//流水线 1 执行访存指令数:
PCR_OPT_DEF(PC7_INST_PIPE1_MEMACC              , 0xc)
//执行向量访存指令数(不含VLENMAS):
PCR_OPT_DEF(PC7_INST_VECTOR_MEMACC             , 0xd)
//向量浮点乘加指令数(含VLENMAS):
PCR_OPT_DEF(PC7_INST_VECTOR_FLOAT_MULS         , 0xe)
//WRCH 请求数:
PCR_OPT_DEF(PC7_REQ_WRCH                       , 0xf)
//全栏栅等待周期数:
PCR_OPT_DEF(PC7_REQ_FENCE_TOTAL_WAIT           , 0x10)
PCR_OPT_DEF(PC7_NULL17                         , 0x11)
PCR_OPT_DEF(PC7_NULL18                         , 0x12)
//从核收到的对 GPR 的 IO 写请求总数:
PCR_OPT_DEF(PC7_REQ_GPR_READ                   , 0x13)
PCR_DEF_END(PC7)

