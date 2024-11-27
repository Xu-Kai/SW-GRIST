#include "lwpf_where.h"
#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif
#include "lwpf_backend.h"
#ifdef __sw_slave__
#include "lwpf_frontend_cpe.h"
#endif
#ifdef __sw_host__
#include "lwpf_frontend_mpe.h"
#endif
