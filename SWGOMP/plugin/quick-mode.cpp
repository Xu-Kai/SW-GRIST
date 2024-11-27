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
#include "swgomp-constants.h"

pass_data pass_data_quickmode = {GIMPLE_PASS,
                               "swgomp-quickmode",
                               OPTGROUP_OMP,
                               TV_NONE,
                               PROP_gimple_any,
                               PROP_gimple_eomp,
                               0,
                               0,
                               0};
class swgomp_quickmode : public gimple_opt_pass {
public:
  swgomp_quickmode(gcc::context *ctx) : gimple_opt_pass(pass_data_quickmode, ctx) {
  }
  static tree my_walk_stmt(gimple_stmt_iterator *gsi, bool *walk_sub, walk_stmt_info *wi){
    auto stmt = gsi ? gsi_stmt(*gsi) : nullptr;
    if (stmt && gimple_code(stmt) == GIMPLE_OMP_TARGET && gimple_omp_target_kind(stmt) == GF_OMP_TARGET_KIND_REGION){
      gomp_target* target_stmt = as_a<gomp_target*>(stmt);
      auto body = gimple_omp_body(stmt);
      bool is_single = true;
      int cnt = 0;
      gomp_parallel *parallel_stmt = nullptr;
      for (auto gsi = gsi_start(body); !gsi_end_p(gsi); gsi_next(&gsi)){
        auto stmt = gsi_stmt(gsi);

        while (gimple_code(stmt) == GIMPLE_BIND) {
          auto gsii = gimple_bind_body((gbind*)stmt);
          if (gimple_seq_first_stmt(gsii) != gimple_seq_last_stmt(gsii)) {
            is_single = false;
          }
          stmt = gimple_seq_first_stmt(gsii);
        }

        if (gimple_code(stmt) != GIMPLE_OMP_RETURN) {
          if (gimple_code(stmt) == GIMPLE_OMP_PARALLEL) {
            parallel_stmt = as_a<gomp_parallel*>(stmt);
          } else {
            is_single = false;
          }
        }
      }
      if (is_single) {
        tree thread_limit = NULL_TREE, num_threads = NULL_TREE, num_teams = NULL_TREE;
        for (tree t = gimple_omp_target_clauses(target_stmt); t; t = OMP_CLAUSE_CHAIN(t)) {
          if (OMP_CLAUSE_CODE(t) == OMP_CLAUSE_THREAD_LIMIT) {
            thread_limit = t;
          } else if (OMP_CLAUSE_CODE(t) == OMP_CLAUSE_NUM_TEAMS) {
            num_teams = t;
          }
        }
        for (tree t = gimple_omp_parallel_clauses(parallel_stmt); t; t = OMP_CLAUSE_CHAIN(t)) {
          if (OMP_CLAUSE_CODE(t) == OMP_CLAUSE_NUM_THREADS) {
            num_threads = t;
          }
        }
        
        if (num_teams == NULL_TREE) {
          num_teams = build_omp_clause(gimple_location(target_stmt), OMP_CLAUSE_NUM_TEAMS);
          OMP_CLAUSE_CHAIN(num_teams) = gimple_omp_target_clauses(target_stmt);
          gimple_omp_target_set_clauses(target_stmt, num_teams);
        }
        OMP_CLAUSE_OPERAND(num_teams, 0) = build_int_cst(intSI_type_node, -1);
        if (thread_limit == NULL_TREE) {
          thread_limit = build_omp_clause(gimple_location(target_stmt), OMP_CLAUSE_THREAD_LIMIT);
          OMP_CLAUSE_CHAIN(thread_limit) = gimple_omp_target_clauses(target_stmt);
          gimple_omp_target_set_clauses(target_stmt, thread_limit);
        }
        OMP_CLAUSE_OPERAND(thread_limit, 0) = num_threads == NULL_TREE ? build_int_cst(intSI_type_node, 0) : OMP_CLAUSE_OPERAND(num_threads, 0);
        if (num_threads == NULL_TREE) {
          num_threads = build_omp_clause(gimple_location(parallel_stmt), OMP_CLAUSE_NUM_THREADS);
          OMP_CLAUSE_CHAIN(num_threads) = gimple_omp_parallel_clauses(parallel_stmt);
          gimple_omp_parallel_set_clauses(parallel_stmt, num_threads);
        }
        OMP_CLAUSE_OPERAND(num_threads, 0) = build_int_cst(intSI_type_node, -1);
        // debug_gimple_stmt(stmt);
        update_stmt(stmt);
      }
    }
    *walk_sub = false;
    return NULL_TREE;
  }

  int unsigned execute(function *f) override {
    basic_block bb;
    const char *fname = IDENTIFIER_POINTER(DECL_ASSEMBLER_NAME(f->decl));
    auto node = cgraph_node::get(f->decl);

    gimple_seq body = gimple_body(f->decl);
    walk_stmt_info info;
    memset(&info, 0, sizeof(info));
    walk_gimple_seq(body, my_walk_stmt, nullptr, &info);

    return 0;
  }
};

void register_quickmode(){
  // register_pass(new swgomp_quickmode(g), PASS_POS_INSERT_AFTER, "omplower", 1);
}