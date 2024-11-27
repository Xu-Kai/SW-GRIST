def calc_bitmap(myid, blev, clev, nthr):
  nl = (nthr+1) >> 1
  nr = nthr - nl
  bm = 0
  
  bm |= 1 << myid
  q = [(myid, clev, nthr)]
  s = 0
  while s < len(q):
    root, clev, nthr = q[s]
    nl = (nthr + 1) >> 1
    nr = nthr - nl
    if nr > 0:
      q.append((root+(1<<5-clev), clev+1, nr))
    if nl > 1:
      q[s] = (root, clev+1, nl)
    else:
      bm |= 1 << root
      s += 1
  return bm
# def spawn(myid, blev, clev, nthr, tid, pc, arg):
#   print(myid, clev, nthr)
#   nl = (nthr + 1) >> 1
#   nr = nthr - nl
#   if nr > 0:
#     spawn(myid+(1<<5-clev), blev, clev+1, nr, tid+1, pc, arg)
#   while nl > 1:
#     clev += 1
#     nthr = nl
#     nl = (nthr+1) >> 1
#     nr = nthr - nl
#     spawn(myid+(1<<5-clev), blev, clev+1, nr, tid+1, pc, arg)
#   #取5-blev到5-clev之间的二进制位取反
# spawn(0, 3, 3, 6, 0, 0, 0)
def print_bm(bm):
  sbm = bin(bm)[::-1]
  for i in range(8):
    print(sbm[i*8:i*8+8])
print_bm(calc_bitmap(0, 0, 0, 60))

# m = [0 for i in range(6)]
# for i in range(64):
#   if i & 1:
#     m[0] |= 1 << i
#   elif i & 2:
#     m[1] |= 1 << i
#   elif i & 4:
#     m[2] |= 1 << i
#   elif i & 8:
#     m[3] |= 1 << i
#   elif i & 16:
#     m[4] |= 1 << i
#   elif i & 32:
#     m[5] |= 1 << i
# for i in m:
#   print(hex(i))