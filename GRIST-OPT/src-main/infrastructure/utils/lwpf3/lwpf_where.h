#ifndef LWPF_WHERE
#define LWPF_WHERE

/* We want to know where we are */
#if defined(__sw_64_thl__)
#ifndef SW5
#define SW5
#endif
#else
#include <crts.h>
/* determine SW7/9 by LDM size */
#if CRTS_MAX_LDM_SIZE == (160 << 10)
#define SW7
#else
#define SW9
#endif /* MAX_LDM_SIZE */
#endif /*__sw_64_thl__*/

#ifdef SW5 /* define the macros for CPEs */

#ifdef SW2
#define __sw_host__
#else
#define __sw_slave__
#endif /* SW2 */

#endif /* SW5 */

#endif /*LWPF_WHERE*/
