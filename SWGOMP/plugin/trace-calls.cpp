
#include <gcc-plugin.h>

#include <langhooks.h>
#include <tree.h>
#include <print-tree.h>
#include <tree-pass.h>
#include <gimple.h>
#include <gimple-pretty-print.h>
#include <gimple-iterator.h>
#include <cgraph.h>
#include <context.h>
#include <stringpool.h>
#include <target.h>
#include <hooks.h>
#include <rtl.h>
#include <genrtl.h>
#include <expr.h>
#include <memmodel.h>
#include <emit-rtl.h>
#include <attribs.h>
#include <string.h>
#include <set>
#include <diagnostic-core.h>
#include <gimple-walk.h>
#include <gimple-ssa.h>
#include "swgomp-constants.h"

pass_data pass_data_swgomp_trace_calls = {
    RTL_PASS, "swgomp-trace-calls", OPTGROUP_NONE, TV_NONE, PROP_rtl, 0, 0};

class swgomp_trace_calls : public rtl_opt_pass {
public:
  swgomp_trace_calls(gcc::context *ctxt)
      : rtl_opt_pass(pass_data_swgomp_trace_calls, ctxt) {}
  virtual bool gate(function *f) { return TARGET_SW_SLAVE; }

  virtual unsigned execute(function *f) {
    for (auto insn = get_insns(); insn; insn = NEXT_INSN(insn)) {
      if (GET_CODE(insn) == CALL_INSN) {
        auto pat = PATTERN(insn);
        rtx pat_call = nullptr;
        if (GET_CODE(pat) == PARALLEL) {
          for (int i = 0; i < XVECLEN(pat, 0); i++) {
            auto para_exp = XVECEXP(pat, 0, i);
            if (GET_CODE(para_exp) == SET) {
              para_exp = XEXP(para_exp, 1);
            }
            if (GET_CODE(para_exp) == CALL) {
              pat_call = para_exp;
            }
          }
        }
        rtx symref = nullptr;
        if (pat_call != nullptr) {
          if (GET_CODE(XEXP(pat_call, 0)) == MEM &&
              GET_CODE(XEXP(XEXP(pat_call, 0), 0)) != SYMBOL_REF) {
            warning_at(INSN_LOCATION(insn), 0,
                       "found non-symbol-ref function call");
          }
        }
      }

      return TODO_df_verify | TODO_df_finish;
    }
    return 0;
  }
};
void register_trace_call(){
  register_pass(new swgomp_trace_calls(g), PASS_POS_INSERT_BEFORE, "sched1",
                1);
}