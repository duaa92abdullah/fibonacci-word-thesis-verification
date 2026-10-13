"""Factor trees of characteristic Sturmian words, Theorem 'factortree' and the table for T_8."""
def charword(ds, length):          # s_{-1}=1, s_0=0, s_j = s_{j-1}^{d_j} s_{j-2}
    a, b = "1", "0"
    for d in ds:
        a, b = b, b * d + a
        if len(b) >= length:
            break
    return b[:length]

def tree_indices(word, n):
    V = {word[i:i + k] for k in range(n + 1) for i in range(len(word) - k + 1)}
    deg = {v: 0 for v in V}
    E = [(v[:-1], v) for v in V if v]
    for u, v in E:
        deg[u] += 1
        deg[v] += 1
    irr = sum(abs(deg[u] - deg[v]) for u, v in E)
    irr2 = sum(abs(deg[u] ** 2 - deg[v] ** 2) for u, v in E)
    M1 = sum(d * d for d in deg.values())
    M2 = sum(deg[u] * deg[v] for u, v in E)
    return len(V), irr, irr2, M1, M2

if __name__ == "__main__":
    n = 8
    print("d, (|V|, irr, irr2, M1, M2) for T_8")
    for d in (1, 2, 3, 4):
        got = tree_indices(charword([d] + [1] * 30, 20000), n)
        formula = ((n + 1) * (n + 2) // 2, 4 * n - 2 * d, 18 * n - 10 * d - 2, 2 * n * n + 8 * n - 4, 2 * n * n + 10 * n - 11 + d)
        print(d, got, "OK" if got == formula else "MISMATCH", formula)
