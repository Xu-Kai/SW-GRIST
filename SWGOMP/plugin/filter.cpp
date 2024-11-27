
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
#include "swgomp-options.hpp"
pass_data pass_data_swgomp_filter = {SIMPLE_IPA_PASS,
                                     "swgomp-filter",
                                     OPTGROUP_NONE,
                                     TV_NONE,
                                     PROP_cfg | PROP_cfglayout,
                                     0,
                                     0,
                                     0,
                                     TODO_dump_symtab};
class swgomp_filter : public simple_ipa_opt_pass {
public:
  swgomp_filter(gcc::context *ctx)
      : simple_ipa_opt_pass(pass_data_swgomp_filter, ctx) {}
  static void mark_declare_target_children(cgraph_node *node) {
    for (cgraph_edge *edge = node->callees; edge; edge = edge->next_callee) {
      auto callee_decl = edge->callee->decl;
      // debug_tree(callee_decl);
      if (edge->callee->definition && callee_decl &&
          !lookup_attribute("omp declare target",
                            DECL_ATTRIBUTES(callee_decl))) {
        DECL_ATTRIBUTES(callee_decl) = make_attribute(
            "omp declare target", "", DECL_ATTRIBUTES(callee_decl));
        inform(DECL_SOURCE_LOCATION(callee_decl),
               "%s is considered for slave code generation\n",
               IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(callee_decl)));
        mark_declare_target_children(edge->callee);
      }
    }
  }
  int unsigned execute(function *) override {
    bool has_slave = false;
    // auto cccmp = [](const char *a, const char *b) { return strcmp(a, b) < 0; };
    // std::set<const char *, decltype(cccmp)> entries(cccmp);
    // if (TARGET_SW_SLAVE) {
    //   // if (getenv("O2ATH_ENTRIES")) {
    //   //   const char *entry_str = xstrdup(getenv("O2ATH_ENTRIES"));
    //   //   int start = 0, iend;
    //   //   for (iend = 0; entry_str[iend]; iend++) {
    //   //     if (entry_str[iend] == ':') {
    //   //       const char *fname = xstrndup(entry_str + start, iend - start);
    //   //       entries.insert(fname);
    //   //       start = iend + 1;
    //   //     }
    //   //   }
    //   //   if (start != iend)
    //   //     entries.insert(xstrndup(entry_str + start, iend - start));
    //   // }
    //   for (auto s : swgomp_entries) {
    //     inform(UNKNOWN_LOCATION,
    //            "%s is considered for generating slave code due to explicit "
    //            "setup\n",
    //            s);
    //   }
    // }
    cgraph_node *node;
    FOR_EACH_FUNCTION(node) {
      auto decl = node->decl;
      const char *fname = IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(decl));
      if (node->offloadable && !lookup_attribute("omp declare target", DECL_ATTRIBUTES(decl)))
        DECL_ATTRIBUTES(decl) =
            make_attribute("omp declare target", "", DECL_ATTRIBUTES(decl));
      if (swgomp_entries.contains(fname)) {
        DECL_ATTRIBUTES(decl) =
            make_attribute("omp declare target", "", DECL_ATTRIBUTES(decl));
        if (TARGET_SW_SLAVE)
          inform(UNKNOWN_LOCATION,
                "%s is considered for generating slave code due to explicit "
                "setup\n",
                fname);
      }
      bool target_entry_p =
          DECL_ATTRIBUTES(decl) &&
          lookup_attribute("omp target entrypoint", DECL_ATTRIBUTES(decl));
      bool target_p =
          DECL_ATTRIBUTES(decl) &&
          lookup_attribute("omp declare target", DECL_ATTRIBUTES(decl));
      if (target_p || target_entry_p) {
        if (node->has_gimple_body_p())
          has_slave = true;
        if (TARGET_SW_SLAVE) {
          if (strncmp("slave_", fname, 6)/* && node->has_gimple_body_p()*/) {
            symtab->change_decl_assembler_name(decl, get_identifier(ACONCAT(("slave_", fname, NULL))));
            // SET_DECL_ASSEMBLER_NAME(
            //     decl, get_identifier(ACONCAT(("slave_", fname, NULL))));
          }
        }
      }
      if (TARGET_SW_HOST) {
        if (target_entry_p) {
          DECL_EXTERNAL(decl) = 1;
          TREE_PUBLIC(decl) = 1;
          node->externally_visible = true;
          symtab->change_decl_assembler_name(decl, get_identifier(ACONCAT(("slave_", fname, NULL))));
          // SET_DECL_ASSEMBLER_NAME(
              // node->decl, get_identifier(ACONCAT(("slave_", fname, NULL))));
          // node->verify_node();
        }
      }
      if (TARGET_SW_SLAVE && (target_p || target_entry_p) &&
          swgomp_recursive_p) {
        mark_declare_target_children(node);
      }
    }
    if (TARGET_SW_SLAVE) {
      FOR_EACH_FUNCTION(node) {
        auto decl = node->decl;
        const char *fname = IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(decl));
        if (swgomp_entries.contains(fname)) {
          DECL_ATTRIBUTES(decl) =
              make_attribute("omp declare target", "", DECL_ATTRIBUTES(decl));
        }
        bool entry =
            lookup_attribute("omp target entrypoint", DECL_ATTRIBUTES(decl));
        bool needed = entry || lookup_attribute("omp declare target",
                                                DECL_ATTRIBUTES(decl));
        if (!node->has_gimple_body_p()) continue;

        if (entry) {
          TREE_PUBLIC(decl) = 1;
          node->externally_visible = true;
        }
        if (needed) {
          node->mark_force_output();
          if (DECL_STATIC_CHAIN(decl)) {
            warning_at(
                DECL_SOURCE_LOCATION(decl), 0,
                "%s is considered as always_inline to avoid static-chain "
                "calls\n",
                fname);
            DECL_DECLARED_INLINE_P(decl) = 1;
            DECL_ATTRIBUTES(decl) = make_attribute(
                "no_instrument_function", "",
                make_attribute("always_inline", "", DECL_ATTRIBUTES(decl)));
            DECL_DISREGARD_INLINE_LIMITS(decl) = 1;
          }
          node->externally_visible = true;
          node->parallelized_function = 0;
          node->resolution = LDPR_RESOLVED_IR;
          // puts(fname);
        } else {
            DECL_EXTERNAL(decl) = 1;
        }
      }
    }
    if (TARGET_SW_HOST && has_slave && swgomp_slave_flag_path) {
      fclose(fopen(swgomp_slave_flag_path, "w"));
    }
    return 0;
  }
};

void register_filter() {
  register_pass(new swgomp_filter(g), PASS_POS_INSERT_BEFORE, "visibility", 1);
}