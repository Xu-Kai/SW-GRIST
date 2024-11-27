import re
import sys
import subprocess
node_re = re.compile(r'("[0-9a-f]*")\s*\[label\s*=\s*"(.*)"\]')
edge_re = re.compile(r'("[0-9a-f]*")\s*->\s*("[0-9a-f]*")\s*\[label\s*=\s*"([0-9a-z]*.[0-9a-z]*%)"\]')
tot_re = re.compile(r'self: ([0-9.]+)% [0-9.]+% [0-9.%]+')
mod_re = re.compile(r'__(\w+)_MOD_(\w+)')
name_re = re.compile(r'^(\w+)')
nodes = []
edges = []
filt = subprocess.Popen("c++filt", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr, universal_newlines=True)
def filt_name(m):
    name = m.group()
    filt.stdin.write(name+"\n")
    filt.stdin.flush()
    name = filt.stdout.readline().strip()
    mm = mod_re.match(name)
    if mm:
        name = "%s::%s" % mm.groups()
    return name
for line in open(sys.argv[1]):
    mn = node_re.fullmatch(line.strip())
    me = edge_re.fullmatch(line.strip())
    if mn:
        name, label = mn.groups()
        tm = tot_re.search(label)
        rat = float(tm.group(1))/100
        rat = (rat ** 0.375)
        label = name_re.sub(filt_name, label)
        fill = "#%02x%02x%02x" % (int((1-rat)*255), int(255), int((1-rat)*255))
        print('%s [label = "%s", fillcolor = "%s", style=filled]' % (name, label, fill))
    elif me:
        s, t, label = me.groups()
        rat = float(label[:-1])/100
        w = rat*7+1
        print('%s -> %s [label = "%s", penwidth=%f]' % (s, t, label, w))
    else:
        print(line.rstrip())
    

