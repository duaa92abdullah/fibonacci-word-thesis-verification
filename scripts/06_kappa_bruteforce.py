"""Theorem 'tw:main': kappa(N,b,j) = max(2, floor(2b/j'))  versus brute force.  (numpy version of the original loop)"""
from math import gcd
import numpy as np

def digit_sums(L, b):
    s = np.zeros(L, dtype=np.int64)
    x = np.arange(L, dtype=np.int64)
    while x.any():
        s += x % b
        x //= b
    return s

def longest_run(mask):
    m = np.concatenate(([0], mask.astype(np.int8), [0]))
    d = np.diff(m)
    st, en = np.where(d == 1)[0], np.where(d == -1)[0]
    return int((en - st).max()) if len(st) else 0

def kappa(N, b, j):
    g = N // gcd(N, b - 1)
    if g == 1:
        return None                   # W periodic: no bound
    jp = (j - 1) % N + 1              # j' in {1..N}, j' = j mod N
    return max(2, (2 * b) // jp)

def brute(N, b, j, L=60000, mmax=200):
    s = digit_sums(L + mmax, b)
    best = 1
    for m in range(1, mmax + 1):
        c = s[m:L + m] - s[:L]
        best = max(best, 1 + longest_run((c % N) == (j % N)) // m)
    return best

if __name__ == "__main__":
    bad = 0
    for N in range(2, 7):
        for b in range(2, 6):
            for j in range(N):
                kp = kappa(N, b, j)
                if kp is None:
                    continue
                kb = brute(N, b, j)
                bad += kb != kp
                print(N, b, j, "theory", kp, "brute", kb, "OK" if kb == kp else "MISMATCH")
    print("mismatches:", bad)
