
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
pass_data pass_data_swgomp_private_name = {GIMPLE_PASS,
                                     "swgomp-private-name",
                                     OPTGROUP_NONE,
                                     TV_NONE,
                                     PROP_cfg | PROP_cfglayout,
                                     0,
                                     0,
                                     0,
                                     0};
class swgomp_private_name : public gimple_opt_pass {
public:
  swgomp_private_name(gcc::context *ctx)
      : gimple_opt_pass(pass_data_swgomp_private_name, ctx) {}
  int unsigned execute(function *f) override {
    // debug_tree(DECL_ASSEMBLER_NAME(f->decl));
    // fprintf(stderr, "%d %d\n", DECL_UID(f->decl), f->funcdef_no);
    const char *oldlabel = IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(f->decl));
    int uid = DECL_UID(f->decl);
    char uid_s[20];
    sprintf(uid_s, "%d", uid);
    if (!strcmp(oldlabel + strlen(oldlabel)-strlen(uid_s), uid_s)) {
      char *newlabel;
      ASM_FORMAT_PRIVATE_NAME(newlabel, IDENTIFIER_POINTER(DECL_NAME(f->decl)), f->funcdef_no);
      // symtab->change_decl_assembler_name(decl, get_identifier(ACONCAT(("slave_", fname, NULL))));
      symtab->change_decl_assembler_name(f->decl, get_identifier(newlabel));
      inform(DECL_SOURCE_LOCATION(f->decl), "found uid labeled function %s, renamed to %s\n", oldlabel, newlabel);
    }
    
    return 0;
  }
};

void register_private_name() {
  register_pass(new swgomp_private_name(g), PASS_POS_INSERT_BEFORE, "omplower", 1);
}