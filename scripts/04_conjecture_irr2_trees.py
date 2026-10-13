"""Conjecture: irr_2(T) < irr_2(S_n) = n(n-1)(n-2) for every tree that is not the star.  Checked for all trees, n=5..12."""
import networkx as nx

for n in range(5, 13):
    star = n * (n - 1) * (n - 2)
    worst = 0
    count = 0
    for T in nx.nonisomorphic_trees(n):
        if max(dict(T.degree()).values()) == n - 1:
            continue
        irr2 = sum(abs(T.degree(u) ** 2 - T.degree(v) ** 2) for u, v in T.edges())
        assert irr2 < star, (n, irr2, star)
        worst, count = max(worst, irr2), count + 1
    print(f"n={n}: {count} non-star trees, max irr_2 = {worst} < {star}")
print("ok")
