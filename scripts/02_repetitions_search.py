"""Search for strong (3, sigma^j)-repetitions in W (Remark after Theorem 'avoidance').
Uses the criterion of Lemma 'arith': block length m, start i  <=>  c_m(n) = s_2(n+m)-s_2(n) = j (mod N)
on 2m consecutive n.  For each m we look at the longest run of such n (fast, numpy)."""
import numpy as np

def digit_sums(L, b=2):
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

def blocks_with_repetition(N, j, mmax=399, L=100000):
    s = digit_sums(L + mmax)
    found = []
    for m in range(1, mmax + 1):
        c = s[m:L + m] - s[:L]
        if longest_run((c % N) == (j % N)) >= 2 * m:     # strong (3, sigma^j)-repetition with block length m
            found.append(m)
    return found

if __name__ == "__main__":
    for N in (3, 4, 5, 6):
        for j in range(N):
            f = blocks_with_repetition(N, j)
            print(f"N={N} j={j}: block lengths m<400 with a (3,sigma^j)-repetition:", f if f else "none")
