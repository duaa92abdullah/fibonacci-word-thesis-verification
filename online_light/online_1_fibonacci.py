# Fibonacci word checks - light version (runs in < 2 s)
import math
def fib_word(n):
    a, b = "1", "0"
    for _ in range(n):
        a, b = b, b + a
    return b
def count(w, u):
    return sum(1 for i in range(len(w) - len(u) + 1) if w[i:i+len(u)] == u)
def fl_alpha(m):                      # floor(m/phi^2), exact
    return (3*m - math.isqrt(5*m*m) - 1)//2 if m else 0
f = fib_word(18)                      # length 4181
# 1. exact count of ones in prefixes
ones = 0; ok = True
for n in range(1, 4000):
    ones += f[n-1] == "1"
    ok &= ones == fl_alpha(n+1)
print("1. |p_n|_1 = floor((n+1)/phi^2), n<4000:", ok)
# 2. the examples of Section 3.4
p = f[:1000]
print("2. ones in p_1000 =", p.count("1"), "(expected 382); 010:", count(p, "010"), "(381); 101:", count(p, "101"), "(145)")
# 3. table of palindrome densities
s4, s5, s6 = fib_word(4), fib_word(5), fib_word(6)
print("3. 00 in s4:", count(s4, "00"), "/", len(s4), "| 101 in s5:", count(s5, "101"), "/", len(s5),
      "| 010 in s5:", count(s5, "010"), "/", len(s5), "| 1001 in s6:", count(s6, "1001"), "/", len(s6))
# 4. s_9 contains all factors of length <= 30, and palindromes per length
def factors(w, L): return {w[i:i+L] for i in range(len(w)-L+1)}
s9 = fib_word(9)
print("4. s_9 (length", len(s9), ") contains all factors of length<=30:", all(factors(s9, L) == factors(f, L) for L in range(1, 31)))
pal = [len({u for u in factors(f, L) if u == u[::-1]}) for L in range(1, 31)]
print("5. palindromes per length:", pal, "->", all(x == (2 if L % 2 else 1) for L, x in zip(range(1, 31), pal)))
print("6. no '11' and no '000':", "11" not in f and "000" not in f)
phi = (1 + 5**0.5)/2
print("7. freq 101 limit 1/phi^4 =", round(1/phi**4, 6), "and 1/phi^3 =", round(1/phi**3, 6), "(the latter is for 001, 100)")
