"""Counts and frequencies of the factors of length 3 in prefixes of the Fibonacci word
(replaces the three duplicated scripts).  Theory: freq(010)=1/phi^2, freq(101)=1/phi^4, freq(001)=freq(100)=1/phi^3."""
import math

def fib_prefix(n):
    a, b = "0", "01"
    while len(b) < n:
        a, b = b, b + a
    return b[:n]

def count_occurrences(word, pattern):
    m = len(pattern)
    return sum(1 for i in range(len(word) - m + 1) if word[i:i + m] == pattern)

PHI = (1 + math.sqrt(5)) / 2
THEORY = {"010": 1 / PHI ** 2, "101": 1 / PHI ** 4, "001": 1 / PHI ** 3, "100": 1 / PHI ** 3}   # NOTE: 101 is 1/phi^4, not 1/phi^3

def analyze(n, pattern):
    w = fib_prefix(n)
    c = count_occurrences(w, pattern)
    th = THEORY[pattern]
    print(f"n={n:6d} {pattern}: count={c:6d}  count/n={c/n:.6f}  count/(n-2)={c/(n-len(pattern)+1):.6f}  theory={th:.6f}  |count/n-theory|={abs(c/n-th):.6f}")

if __name__ == "__main__":
    assert abs(sum(THEORY.values()) - 1) < 1e-12          # the four length-3 factors have total frequency 1
    for n in (100, 200, 500, 1000, 2000, 20000):
        for pat in ("010", "101"):
            analyze(n, pat)
