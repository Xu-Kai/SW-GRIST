
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
#include <gimplify.h>
#include <stdio.h>
#include "swgomp-options.hpp"
#define GFC_DESCRIPTOR_TYPE_P(node) TYPE_LANG_FLAG_1(node)
pass_data pass_data_swgomp_protect_optional = {GIMPLE_PASS,
                                     "swgomp-protect-optional",
                                     OPTGROUP_NONE,
                                     TV_NONE,
                                     PROP_gimple_any,
                                     0,
                                     0,
                                     0,
                                     0};
class swgomp_protect_optional : public gimple_opt_pass {
public:
  swgomp_protect_optional(gcc::context *ctx)
      : gimple_opt_pass(pass_data_swgomp_protect_optional, ctx) {}
  bool gate(function *f) override {
    return lang_GNU_Fortran();
  }
  int unsigned execute(function *f) override {
      walk_stmt_info wi;
      push_gimplify_context();
      walk_gimple_seq(f->gimple_body, protect_optional, NULL, &wi);
      pop_gimplify_context(f->gimple_body);
    return 0;
  }
  // walk_tree_fn f;
  static tree protect_optional(gimple_stmt_iterator *gsi, bool *walk_sub, walk_stmt_info *wi){
    gimple *stmt = gsi_stmt(*gsi);
    if (gimple_code(stmt) == GIMPLE_ASSIGN && gimple_assign_rhs_code(stmt) == COMPONENT_REF) {
      tree rhs = gimple_assign_rhs1(stmt);
      tree decl = walk_tree(&rhs, find_decl, NULL, NULL);
      /* 改动：用来判断一个传入的参数是不是optional的方法是基于输入是指针还是引用猜测*/
      /*POINTER_TYPE_P(t)对REFERENCE_TYPE和POINTER_TYPE均有效，故只好基于TREE_CODE判断*/
      if (decl && TREE_CODE(TREE_TYPE(decl)) == POINTER_TYPE) { 
        tree desc_type = TREE_TYPE(TREE_TYPE(decl));
        if (GFC_DESCRIPTOR_TYPE_P(desc_type)) {
          auto rtype = copy_node(TREE_TYPE(rhs));
          if (POINTER_TYPE_P(rtype) || INTEGRAL_TYPE_P(rtype) || FLOAT_TYPE_P(rtype)) {
            auto nonnull = build2(NE_EXPR, boolean_type_node, decl, null_pointer_node);
            auto cond = build3(COND_EXPR, rtype, nonnull, rhs, build1(CONVERT_EXPR, rtype, build_int_cst(intSI_type_node, 0)));
            gimple_seq st = NULL;
            gimplify_assign(gimple_assign_lhs(stmt), cond, &st);
            tree block = make_node(BLOCK);
            gbind *bind = gimple_build_bind(NULL, NULL, block);
            gsi_replace_with_seq(gsi, st, false);
          }
          *walk_sub = true;
        }
      }
    }
    *walk_sub = false;
    return NULL_TREE;
  }
  static tree find_decl(tree *t, int *walk_sub, void *) {
    if (TREE_CODE(*t) == COMPONENT_REF 
        && TREE_CODE(TREE_OPERAND(*t, 0)) == MEM_REF 
        && TREE_CODE(TREE_OPERAND(TREE_OPERAND(*t, 0), 0)) == PARM_DECL) {
      walk_sub = 0;
      return TREE_OPERAND(TREE_OPERAND(*t, 0), 0);
    }
    return NULL_TREE;
  }
};

void register_protect_optional() {
  register_pass(new swgomp_protect_optional(g), PASS_POS_INSERT_AFTER, "omplower", 1);
}
