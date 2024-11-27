
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
#include <diagnostic-core.h>
#include <gimple-walk.h>
#include <gimple-ssa.h>
#include "swgomp-options.hpp"
pass_data pass_data_swgomp_slave_mangle = {RTL_PASS,
                                     "swgomp-slave-mangle",
                                     OPTGROUP_NONE,
                                     TV_NONE,
                                     PROP_cfg | PROP_cfglayout,
                                     0,
                                     0,
                                     0,
                                     TODO_verify_all};
class swgomp_slave_mangle : public rtl_opt_pass {
public:
  swgomp_slave_mangle(gcc::context *ctx)
      : rtl_opt_pass(pass_data_swgomp_slave_mangle, ctx) {}
  int unsigned execute(function *f) override {
    const char *fname = IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(f->decl));
    // cgraph_node::get(f->decl)->verify_node();
    if (TARGET_SW_SLAVE){
      if (strncmp("slave_", fname, 6)) {
        SET_DECL_ASSEMBLER_NAME(
            f->decl, get_identifier(ACONCAT(("slave_", fname, NULL))));
      }
    }
    cgraph_node::get(f->decl)->analyze();
    return 0;
  }
};

void register_slave_mangle() {
  // register_pass(new swgomp_slave_mangle(g), PASS_POS_INSERT_AFTER, "dfinish", 1);
}