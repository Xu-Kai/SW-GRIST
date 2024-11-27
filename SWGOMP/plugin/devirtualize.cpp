
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
pass_data pass_data_swgomp_devirtualize = {GIMPLE_PASS,
                                           "swgomp-devirtualize-fortran",
                                           OPTGROUP_OMP,
                                           TV_NONE,
                                           PROP_gimple_any,
                                           PROP_gimple_eomp,
                                           0,
                                           0,
                                           0};
class swgomp_devirtualize : public gimple_opt_pass {
public:
  swgomp_devirtualize(gcc::context *ctx)
      : gimple_opt_pass(pass_data_swgomp_devirtualize, ctx) {}

  virtual bool gate(function *fun) override { return TARGET_SW_SLAVE; }
  int unsigned execute(function *f) override {
    // basic_block bb;
    // const char *fname = IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(f->decl));
    // auto node = cgraph_node::get(f->decl);
    for (auto gsi = gsi_start(f->gimple_body); !gsi_end_p(gsi);
         gsi_next(&gsi)) {
      auto stmt = gsi_stmt(gsi);
      if (stmt && gimple_code(stmt) == GIMPLE_ASSIGN &&
          gimple_assign_rhs_code(stmt) == COMPONENT_REF) {
        auto rhs = gimple_assign_rhs1(stmt);
        auto field = TREE_OPERAND(rhs, 1);
        auto object = TREE_OPERAND(rhs, 0);
        if (DECL_NAME(field) &&
            !strcmp(IDENTIFIER_POINTER(DECL_NAME(field)), "_vptr")) {
          auto vtype = TREE_TYPE(TREE_TYPE(field));
          const char *name = IDENTIFIER_POINTER(TYPE_NAME(vtype));
          const char *ns;
          if (TYPE_CONTEXT(vtype)) {
            ns = ACONCAT(("__",
                          IDENTIFIER_POINTER(DECL_NAME(TYPE_CONTEXT(vtype))),
                          "_MOD_", NULL));
          } else {
            ns = "";
          }
          const char *vtab_asm_name = ACONCAT((ns, "__vtab", name + 7, NULL));
          puts(vtab_asm_name);
          varpool_node *vn =
              varpool_node::get_for_asmname(get_identifier(vtab_asm_name));
          if (!vn) {
            auto vtab_decl = build_decl(UNKNOWN_LOCATION, VAR_DECL,
                                        get_identifier(vtab_asm_name), vtype);
            SET_DECL_ASSEMBLER_NAME(vtab_decl, get_identifier(vtab_asm_name));
            vn = varpool_node::get_create(vtab_decl);
          }
          // debug_tree(TYPE_CONTEXT(vtype));
          gimple_assign_set_rhs_code(stmt, ADDR_EXPR);
          gimple_assign_set_rhs1(
              stmt, build_fold_addr_expr_loc(gimple_location(stmt), vn->decl));
          warning_at(gimple_location(stmt), 0,
                     "aggressively devirtualized _vptr for %s::%s to %s", name,
                     ns, vtab_asm_name);
        }
      }
    }
    // }
    return 0;
  }
};
void register_devirtualize() {
  register_pass(new swgomp_devirtualize(g), PASS_POS_INSERT_BEFORE, "ssa", 1);
}