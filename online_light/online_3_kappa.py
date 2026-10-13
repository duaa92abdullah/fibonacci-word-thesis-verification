# Theorem tw:main: kappa(N,b,j) = max(2, floor(2b/j'))  versus brute force (small cases) - light version
import numpy as np
from math import gcd
def digit_sums(L, b):
    s = np.zeros(L, dtype=np.int64); x = np.arange(L, dtype=np.int64)
    while x.any():
        s += x % b; x //= b
    return s
def longest_run(mask):
    m = np.concatenate(([0], mask.astype(np.int8), [0])); d = np.diff(m)
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    return int((en-st).max()) if len(st) else 0
bad = 0
for N, b in ((2, 2), (3, 2), (3, 3), (4, 3), (5, 2)):
    g = N // gcd(N, b-1)
    L = 4000; mmax = 40
    s = digit_sums(L + mmax, b)
    for j in range(N):
        jp = (j-1) % N + 1
        theory = max(2, (2*b)//jp)
        best = 1
        for m in range(1, mmax+1):
            c = s[m:L+m] - s[:L]
            best = max(best, 1 + longest_run((c % N) == (j % N))//m)
        bad += best != theory
        print("N=%d b=%d j=%d theory=%d brute=%d %s" % (N, b, j, theory, best, "OK" if best == theory else "MISMATCH"))
print("mismatches:", bad)
