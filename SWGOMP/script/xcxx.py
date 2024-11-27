#!/usr/bin/env python3
import sys
import os
import subprocess
import tempfile
import re
import shutil
import time
c = False
v = False
k = False
vd = False
o = None
src = None
modulepath = None
compiler_set = {
    "xfort": ("sw9gfortran", "mpif90"),
    "xcxx" : ("sw9g++", "mpicxx"),
    "xcc"  : ("sw9gcc", "mpicc"),
}
script_name = os.path.basename(sys.argv[0])
if script_name.endswith(".py"):
    script_name = script_name[:-3]
with_mpi = script_name.startswith("mpi")
if with_mpi:
    script_name = script_name[3:]
with_fortran = script_name == "xfort"
serialc, mpic = compiler_set[script_name]
hostc = mpic if with_mpi else serialc
slavec = serialc
link = mpic if with_mpi else serialc
args_link = ["-Wl,--wrap,GOMP_target_ext,--wrap,athread_init,--wrap,__real_athread_spawn,--wrap,athread_join"]
args = []
prefix = os.path.dirname(os.path.abspath(sys.argv[0]))
plugin = os.path.join(prefix, "..", "plugin", "swgomp.so")
plugin_arg = "-fplugin-arg-swgomp-"
lib = os.path.join(prefix, "..", "lib", "libswgomp.a")
SWGOMP_RE_F90 = re.compile(r"\s*\!\$\s*swgomp\s+entry\((.*)\)\s*")
SWSLAVE_RE = re.compile(r"\s*#ifdef\s+__sw_slave__")
i = 1
argv = list(filter(lambda x: not x.startswith("-Wslave,") and not x.startswith("-Whost",), sys.argv))
args_slave = []
args_host = []
while i < len(argv):
    arg = argv[i].strip()
    #print(i, arg)
    if arg == "-c":
        c = True
    elif arg == "-v":
        v = True
    elif arg == "-vv":
        args.append("-v")
        v = True
    elif arg == "-k":
        k = True
    elif arg == "-o":
        o = argv[i+1]
        i += 1
    elif arg == "-vd":
        vd = True
    elif arg == "-fno-swgomp-recursive":
        args.append("%sno-recursive" % plugin_arg)
    elif arg == "-save-temps":
        args.append("-save-temps=obj")
        k = True
    elif arg.startswith("-fdump-"):
        k = True
        args.append(arg)
    elif arg == "-I.":
        args.append("-I%s" % os.path.abspath('.'))
    elif arg.startswith("-J"):
        modulepath = arg[2:]
    elif os.path.splitext(arg)[-1] in [".f90", ".F90", ".F", ".f", ".in", ".c", ".cpp", ".C", ".cc"]:
        src = arg
        args.append(arg)
    else:
        args.append(arg)
    i += 1
for slf in filter(lambda x: x.startswith("-Wslave,"), sys.argv):
    args_slave += slf.split(",")[1:]
for hof in filter(lambda x: x.startswith("-Whost,"), sys.argv):
    args_host += hof.split(",")[1:]
def run(args, fatal=True):
    if v:
        vargs = [x if len(x.split()) == 1 else "'%s'" % x.replace("'", "\\'") for x in args]
        print(" ".join(vargs), file=sys.stderr)
    ret = subprocess.run(args, stderr=sys.stderr, stdout=sys.stdout).returncode
    if fatal and ret:
        vargs = [x if len(x.split()) == 1 else "'%s'" % x.replace("'", "\\'") for x in args]
        print(" ".join(vargs) + " returns non-zero status %d" % ret, file=sys.stderr)
        sys.exit(ret)
    return ret
def check(args):
    if v:
        print(args, file=sys.stderr)
    return subprocess.check_output(args, stderr=sys.stderr).decode()
if o is None:
    if src is not None:
        if c:
            o = os.path.splitext(os.path.basename(src))[0] + ".o"
        else:
            o = os.path.join(os.getcwd(), "a.out")
if modulepath is None:
    modulepath = "."
if src is None:
    sys.exit(run([link, "-mhybrid", "-fopenmp"] + args_link + argv[1:] + [lib]))
win = False
start_time = time.time()
try:
    if not k:
        work = tempfile.mkdtemp(prefix="swgomp-")
    else:
        work = o + ".dir"
        if os.path.exists(work):
            shutil.rmtree(work)
        os.mkdir(work)
    args_host += ["-dumpdir", os.path.join(work, ""), "-dumpbase", os.path.basename(src) + ".host"]
    args_slave += ["-dumpdir", os.path.join(work, ""), "-dumpbase", os.path.basename(src) + ".slave"]
    if with_fortran:
        args_host.append("-J%s" % work)
        args.append("-I%s" % os.path.join(prefix, "..", "include/fortran/"))
    else:
        args_slave.append("-I%s" % os.path.join(prefix, "..", "include/slave/"))
    if with_mpi:
        args_slave.append("-I/usr/sw/mpi/mpi_current/include/")
    obase = os.path.basename(os.path.splitext(o)[0])
    host_o = os.path.join(work, obase + ".host.o")
    slave_o = os.path.join(work, obase + ".slave.o")
    slave_flag = os.path.join(work, "slave_flag")
    args_host.append("%sslave-flag=%s" % (plugin_arg, slave_flag))
    os.environ["SLAVE_FLAG"] = slave_flag
    entries = []
    force_slave = False
    if src is not None:
        for line in open(src):
            if SWGOMP_RE_F90.fullmatch(line):
                entries = entries + list(map(str.strip, SWGOMP_RE_F90.fullmatch(line).group(1).split(",")))
            if SWSLAVE_RE.match(line):
                force_slave = True
    os.environ["O2ATH_ENTRIES"] = ":".join(entries)
    args += ['%sentry=%s' % (plugin_arg, ent) for ent in entries]
    if entries:
        print("making the plugin to consider %s for slave code generation" % str(entries), file=sys.stderr)
    run([hostc, "-mhost", "-c", "-fplugin=%s" % plugin, "-fopenmp", "-D__sw_host__"] + args + args_host + ["-o", host_o])

    if os.path.exists(slave_flag) or entries or force_slave:
        run([slavec, "-mslave", "-c", "-fplugin=%s" % plugin, "-fopenmp", "-D__sw_slave__"] + args + args_slave + ["-o", slave_o])
        if c:
            run(["swld", "-r", "-z", "muldefs", host_o, slave_o, "-o", o])
        else:
            run([hostc, "-mhybrid", "-fopenmp", host_o, slave_o, "-o", o] + args_link + [lib])
    else:
        if c:
            run(["mv", host_o, o])
        else:
            run([hostc, "-mhybrid", "-fopenmp", host_o, "-o", o, lib] + args_link)
    if with_fortran:
        for f in os.listdir(work):
            if f.endswith(".mod"):
                if not os.path.exists(os.path.join(modulepath, f)) or open(os.path.join(work, f), "rb").read() != open(os.path.join(modulepath, f), "rb").read():
                    ov = v
                    v = True
                    run(["cp", "-f", os.path.join(work, f), os.path.join(modulepath, f)])
                    v = ov
    win = True
finally:
    if not win:
        if os.path.exists(src + ".o") and os.stat(src + ".o").st_mtime > start_time:
            os.remove(src + ".o")
        if os.path.exists(o):
            os.remove(o)
    if not k:
        shutil.rmtree(work)
    else:
        print("temp files saved in %s" % work)
