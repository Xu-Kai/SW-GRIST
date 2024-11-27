
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
// #include "gfortran.h"
#include <diagnostic-core.h>
#include <gimple-walk.h>
#include <gimple-ssa.h>
#include <c-family/c-pragma.h>
#include "swgomp-constants.h"
#include "swgomp-options.hpp"
static attribute_spec attribute_sw_64_host = {"sw_64_host", 0,     0,    false,
                                              false,        false, NULL, false};
static attribute_spec attribute_sw_64_slave = {
    "sw_64_slave", 0, 0, false, false, false, NULL, false};

static void initialize_attributes(void *, void *) {
  register_attribute(&attribute_sw_64_host);
  register_attribute(&attribute_sw_64_slave);
}

enum cpp_ttype (*xpragma_lex) (tree *, location_t *loc) = nullptr;
static void handle_swgomp_entry(cpp_reader *ARG_UNUSED(dummy)){
  enum cpp_ttype token;
  location_t loc;
  tree x;
  token = xpragma_lex(&x, &loc);
  if (token == CPP_NAME && !strcmp(IDENTIFIER_POINTER(x), "entries")) {
    token = xpragma_lex(&x, &loc);
    if (token == CPP_OPEN_PAREN) {
      
      do {
        token = xpragma_lex(&x, &loc);
        if (token != CPP_NAME) {
          error_at(loc, "invalid swgomp entries pragma, function names needed here.");
        }
        swgomp_entries.add(xstrdup(IDENTIFIER_POINTER(x)));
        token = xpragma_lex(&x, &loc);
      } while (token == CPP_COMMA);
      if (token != CPP_CLOSE_PAREN) {
        error_at(loc, "invalid swgomp entries pragma, \")\" or \",\" needed here.");
      }
      token = xpragma_lex(&x, &loc);
      if (token != CPP_EOF) {
        warning_at(loc, 0, "extra characters in swgomp pragma ignored");
      }
    } else {
      warning_at(loc, 0, "straying swgomp entries pragma.");
    }

  } else {
    if (token != CPP_EOF) {
      warning_at(loc, 0, "invalid swgomp pragma ignored, only support {entries}.");
    }
  }

}

static void register_pragmas(void *ev, void *dat){
  xpragma_lex = (enum cpp_ttype (*) (tree *, location_t *loc))dlsym(RTLD_NEXT, "_Z10pragma_lexPP9tree_nodePj");
  void (*register_pragma) (const char *space, const char *name,
                               pragma_handler_1arg handler) = (void (*) (const char *space, const char *name,
                               pragma_handler_1arg handler))dlsym(RTLD_NEXT, "_Z17c_register_pragmaPKcS0_PFvP10cpp_readerEq");
  if (xpragma_lex && register_pragma) register_pragma(NULL, "swgomp", handle_swgomp_entry);
}
decltype(targetm) targetm_orig;

void dump_decls(void *gcc_data, void *user_data) {
  varpool_node *node;
  FOR_EACH_DEFINED_VARIABLE(node) { node->dump(stdout); }
}
int plugin_is_GPL_compatible = 1;
extern void register_quickmode();
extern void register_filter();
extern void register_devirtualize();
extern void register_trace_call();
extern void register_private_name();
extern void register_slave_mangle();
extern void register_swgomp_lower();
extern void register_protect_optional();
bool swgomp_recursive_p = true;
hash_set<const char*, nofree_string_hash> swgomp_entries;
const char *swgomp_slave_flag_path = nullptr;
int plugin_init(struct plugin_name_args *plugin_info,
                struct plugin_gcc_version *version) {
  targetm_orig = targetm;
  for (int i = 0; i < plugin_info->argc; i ++) {
    const char *k = plugin_info->argv[i].key;
    const char *v = plugin_info->argv[i].value;
    if (!strcmp(k, "slave-flag"))
      swgomp_slave_flag_path = xstrdup(v);
    else if (!strcmp(k, "recursive"))
      swgomp_recursive_p = true;
    else if (!strcmp(k, "no-recursive"))
      swgomp_recursive_p = false;
    else if (!strcmp(k, "entry"))
      swgomp_entries.add(xstrdup(v));
  }
  register_callback("swgomp", PLUGIN_ATTRIBUTES, initialize_attributes, nullptr);
  register_callback("swgomp", PLUGIN_PRAGMAS, register_pragmas, nullptr);
  register_private_name();
  register_devirtualize();
  register_quickmode();
  register_filter();
  register_slave_mangle();
  register_trace_call();
  register_swgomp_lower();
  register_protect_optional();
  return 0;
}