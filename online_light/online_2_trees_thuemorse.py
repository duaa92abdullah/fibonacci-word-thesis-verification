# Factor trees (Theorem factortree, table T_8) and complexity of W for N=3 - light version
def charword(ds, length):
    a, b = "1", "0"
    for d in ds:
        a, b = b, b*d + a
        if len(b) >= length: break
    return b[:length]
def tree_indices(word, n):
    V = {word[i:i+k] for k in range(n+1) for i in range(len(word)-k+1)}
    deg = {v: 0 for v in V}
    E = [(v[:-1], v) for v in V if v]
    for u, v in E:
        deg[u] += 1; deg[v] += 1
    irr = sum(abs(deg[u]-deg[v]) for u, v in E)
    irr2 = sum(abs(deg[u]**2-deg[v]**2) for u, v in E)
    return len(V), irr, irr2, sum(d*d for d in deg.values()), sum(deg[u]*deg[v] for u, v in E)
n = 8
for d in (1, 2, 3, 4):
    got = tree_indices(charword([d] + [1]*30, 3000), n)
    formula = ((n+1)*(n+2)//2, 4*n-2*d, 18*n-10*d-2, 2*n*n+8*n-4, 2*n*n+10*n-11+d)
    print("d =", d, got, "OK" if got == formula else "MISMATCH")
# generalized Thue-Morse, N=3, b=2
def sb(n):
    s = 0
    while n: s += n & 1; n >>= 1
    return s
N = 3
w = [sb(i) % N for i in range(2**13)]
pk = [len({tuple(w[i:i+k]) for i in range(len(w)-k+1)}) for k in range(1, 41)]
print("p(40) =", pk[39], "(expected 303);  differences:", sorted({pk[i+1]-pk[i] for i in range(39)}), "(expected 6 and 9)")
