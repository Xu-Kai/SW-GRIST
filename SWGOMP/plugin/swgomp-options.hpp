#include <gcc-plugin.h>
extern bool swgomp_recursive_p;
extern hash_set<const char*, nofree_string_hash> swgomp_entries;
extern const char *swgomp_slave_flag_path;