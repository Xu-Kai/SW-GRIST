module omnicopy_mod
  interface omnicopy
    module procedure omnicopy_real8
    module procedure omnicopy_real4
    module procedure omnicopy_integer4
    module procedure omnicopy_integer8
    module procedure omnicopy_s_real8
    module procedure omnicopy_s_real4
    module procedure omnicopy_s_integer4
    module procedure omnicopy_s_integer8
  end interface
contains
#define EXP(x) x
#define DTYPE real
#define KIND 8
#include "omnicopy.inc"
#undef DTYPE
#undef KIND
#define DTYPE integer
#define KIND 8
#include "omnicopy.inc"
#undef DTYPE
#undef KIND
#define DTYPE real
#define KIND 4
#include "omnicopy.inc"
#undef DTYPE
#undef KIND
#define DTYPE integer
#define KIND 4
#include "omnicopy.inc"
#undef DTYPE
#undef KIND
end module omnicopy_mod