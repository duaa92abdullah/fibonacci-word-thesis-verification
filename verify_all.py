#!/usr/bin/env python3
"""
Computational verification of the claims in the thesis
"Combinatorics on Words and Symbolic Dynamics: ... Fibonacci Word Density ..."
(D. Abdullah, MIPT, 2026).

Run:   python verify_all.py            (all sections)
       python verify_all.py A C        (only sections A and C)

Every check prints PASS / FAIL.  Exit code is 1 if any check fails.
Sections
  A  Fibonacci word: exact counts, discrepancy, balance, local structure (Ch.3, 3.4)
  B  Examples/tables of Ch.3 (densities of factors, palindromes, letter counts)
  C  Beatty description (Lemma mechanical)
  D  Framed words q_m and error formula (Thm thm1Asymptotic)
  E  Millin-type telescoping series and Fibonacci identities (3.5)
  F  Palindromic prefixes: Fibonacci, ternary word y, characteristic Sturmian words,
     m-bonacci words (3.6)
  G  Generalized Thue-Morse words: strong (k, sigma^j)-repetitions, kappa(N,b,j),
     complexity observation (3.8)
  H  Factor trees of Sturmian words and irregularity indices (Ch.5)
  I  Trees: Prop. pron1, conjecture irr_2 < irr_2(star), Lemma lemn01, A(T)>=leaves
  J  Two-dimensional recursion T(u,v) (Prop. newprove)
  K  Concatenation of base-b expansions of n! (Ch.2, finite sanity check only)
  L  Statistics of Ch.4 recomputed from the SUMMARY numbers in the thesis tables
"""
import sys, math, random
from fractions import Fraction
from decimal import Decimal, getcontext

import numpy as np
import networkx as nx

getcontext().prec = 80
sys.set_int_max_str_digits(0)
RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append((name, bool(cond)))
    print(("PASS  " if cond else "FAIL  ") + name + ((" | " + detail) if detail else ""))


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def fib_list(n):
    F = [0, 1]
    while len(F) <= n:
        F.append(F[-1] + F[-2])
    return F


F = fib_list(1200)
PHI = (1 + Decimal(5).sqrt()) / 2
PSI = (1 - Decimal(5).sqrt()) / 2


def s_word(n):
    """s_n = Theta^n(0), s_{-1}=1, s_0=0, s_n = s_{n-1} s_{n-2}."""
    a, b = "1", "0"
    for _ in range(n):
        a, b = b, b + a
    return b


def theta(w):
    return "".join("01" if c == "0" else "0" for c in w)


def fl_alpha(m):
    """floor(m * alpha), alpha = 1/phi^2 = (3 - sqrt5)/2, exact integer arithmetic."""
    if m == 0:
        return 0
    return (3 * m - math.isqrt(5 * m * m) - 1) // 2


def fl_phi2(m):
    """floor(m * phi^2) = floor(m (3+sqrt5)/2)."""
    return (3 * m + math.isqrt(5 * m * m)) // 2


def fl_phi(m):
    """floor(m * phi) = floor(m (1+sqrt5)/2)."""
    return (m + math.isqrt(5 * m * m)) // 2


def count_overlapping(w, u):
    return sum(1 for i in range(len(w) - len(u) + 1) if w[i:i + len(u)] == u)


def pal_prefix_lengths(w, maxL):
    out = []
    for L in range(1, maxL + 1):
        if w[0] == w[L - 1]:
            p = w[:L]
            if p == p[::-1]:
                out.append(L)
    return out


def pal_factors(w):
    S = set()
    n = len(w)
    for c in range(n):
        l = r = c
        while l >= 0 and r < n and w[l] == w[r]:
            S.add(w[l:r + 1]); l -= 1; r += 1
        l, r = c, c + 1
        while l >= 0 and r < n and w[l] == w[r]:
            S.add(w[l:r + 1]); l -= 1; r += 1
    return S


# ----------------------------------------------------------------------------
# A  Fibonacci word
# ----------------------------------------------------------------------------
def sec_A():
    print("\n== A. Fibonacci word: Theorem mainre, Prop. balance, Lemma local ==")
    # s_n = Theta^n(0)
    w = "0"
    ok = True
    for n in range(0, 14):
        ok &= (w == s_word(n))
        w = theta(w)
    check("A1 s_n = Theta^n(0), n<=13", ok)
    ok = all(len(s_word(n)) == F[n + 2] and s_word(n).count("0") == F[n + 1]
             and s_word(n).count("1") == F[n] for n in range(0, 25))
    check("A2 |s_n|=F_{n+2}, |s_n|_0=F_{n+1}, |s_n|_1=F_n (n<=24)", ok)

    f = s_word(23)  # length F_25 = 75025
    N = 30000
    alpha = (3 - Decimal(5).sqrt()) / 2
    ones = 0
    exact = True
    lo, hi = Decimal(10), Decimal(-10)
    for n in range(1, N + 1):
        ones += f[n - 1] == "1"
        if ones != fl_alpha(n + 1):
            exact = False
        d = Decimal(ones) - n * alpha
        lo, hi = min(lo, d), max(hi, d)
    check("A3 |p_n|_1 = floor((n+1)/phi^2) for n<=30000", exact)
    inv_phi = 1 / PHI
    check("A4 -1/phi < |p_n|_1 - n/phi^2 < 1/phi^2 (n<=30000)",
          lo > -inv_phi and hi < alpha, f"inf={float(lo):.6f}, sup={float(hi):.6f}")
    check("A5 constants sharp: inf -> -1/phi, sup -> 1/phi^2",
          lo < -inv_phi + Decimal("0.001") and hi > alpha - Decimal("0.001"))

    # mechanical word
    c = "".join(str(fl_alpha(n + 1) - fl_alpha(n)) for n in range(1, 5001))
    check("A6 f[n]=floor((n+1)a)-floor(na) (5000 letters)", c == f[:5000])

    # balance of all factors of length m
    pre = [0]
    for ch in f[:8000]:
        pre.append(pre[-1] + (ch == "1"))
    ok = True
    for m in range(1, 120):
        lo_m = fl_alpha(m)
        for i in range(0, 8000 - m):
            k = pre[i + m] - pre[i]
            if k not in (lo_m, lo_m + 1):
                ok = False
                break
        if not ok:
            break
    check("A7 Prop. balance: every factor of length m<120 has floor/ceil(m/phi^2) ones", ok)
    zeros_ok = all(((m - fl_alpha(m)) in (m - fl_alpha(m), m - fl_alpha(m) - 1)) for m in range(1, 50))
    # factor of length 29: 17 or 18 zeros
    zs = {sum(1 for ch in f[i:i + 29] if ch == "0") for i in range(0, 5000)}
    check("A8 factors of length 29 have 17 or 18 zeros", zs == {17, 18}, str(sorted(zs)))

    # local structure
    check("A9 '11' and '000' do not occur", "11" not in f and "000" not in f)
    f3 = {f[i:i + 3] for i in range(len(f) - 2)}
    check("A10 factors of length 3 = {001,010,100,101}", f3 == {"001", "010", "100", "101"})
    # complexity p(n)=n+1
    ok = all(len({f[i:i + n] for i in range(0, 20000)}) == n + 1 for n in range(1, 41))
    check("A11 factor complexity p(n)=n+1, n<=40", ok)


# ----------------------------------------------------------------------------
# B  examples and tables of Chapter 3
# ----------------------------------------------------------------------------
def sec_B():
    print("\n== B. Examples and tables of Section 3.4 ==")
    f = s_word(23)
    p1000 = f[:1000]
    c010, c101 = count_overlapping(p1000, "010"), count_overlapping(p1000, "101")
    check("B1 occurrences of 010 in p_1000 equal 381 (thesis)", c010 == 381, f"computed {c010}")
    check("B2 occurrences of 101 in p_1000 equal 145 (thesis)", c101 == 145, f"computed {c101}")
    check("B3 ones in p_1000 = floor(1001/phi^2) = 382 (thesis counts 381 occurrences of 010; "
          "the 382nd '1' is at the boundary)", p1000.count("1") == fl_alpha(1001) == 382,
          f"ones={p1000.count('1')}, last letters={p1000[-3:]}, first letters={p1000[:3]}")
    # limits of the factor frequencies on a long prefix
    big = s_word(30)  # length F_32 = 2178309
    n = len(big)
    f010 = count_overlapping(big[:300000], "010") / 300000
    f101 = count_overlapping(big[:300000], "101") / 300000
    check("B4 freq(010) -> 1/phi^2", abs(f010 - 0.381966) < 1e-4, f"{f010:.6f}")
    check("B5 freq(101) -> 1/phi^4 = 0.145898", abs(f101 - 0.145898) < 1e-4, f"{f101:.6f}")
    # table Example thm6.2
    s4, s5, s6 = s_word(4), s_word(5), s_word(6)
    check("B6 D(00,s_4)=1/8", Fraction(count_overlapping(s4, "00"), len(s4)) == Fraction(1, 8))
    check("B7 D(101,s_5)=1/13 and D(010,s_5)=4/13",
          Fraction(count_overlapping(s5, "101"), len(s5)) == Fraction(1, 13)
          and Fraction(count_overlapping(s5, "010"), len(s5)) == Fraction(4, 13))
    check("B8 D(1001,s_6)=4/21", Fraction(count_overlapping(s6, "1001"), len(s6)) == Fraction(4, 21))
    # Example: factors of length 3 of s_4
    s4w = s_word(4)
    check("B9 s_4=01001010 and F_3(s_4)={001,010,100,101}",
          s4w == "01001010" and {s4w[i:i + 3] for i in range(6)} == {"001", "010", "100", "101"})
    # palindromic complexity table
    s9 = s_word(9)
    long_f = s_word(23)[:20000]
    ok = True
    for L in range(1, 31):
        P9 = {x for x in pal_factors(s9) if len(x) == L}
        Pl = {x for x in pal_factors(long_f[:3000]) if len(x) == L}
        expected = 1 if L % 2 == 0 else 2
        ok &= (len(P9) == expected == len(Pl))
    check("B10 palindromic complexity 1 (even) / 2 (odd), L<=30, s_9 contains all", ok)
    allf = {long_f[i:i + 30] for i in range(20000 - 29)}
    minlen = next(n for n in range(30, 400) if allf <= {long_f[i:i + 30] for i in range(n - 29)})
    check("B11a s_9 (length 89) contains all FACTORS of length 30 (statement in thesis)",
          allf <= {s9[i:i + 30] for i in range(len(s9) - 29)}, f"shortest prefix containing all 31 factors of length 30 has length {minlen}")
    s10 = s_word(10)
    check("B11b s_10 (length 144) contains all factors of length 30", allf <= {s10[i:i + 30] for i in range(len(s10) - 29)})
    ok = all({x for x in pal_factors(s9) if len(x) == L} == {x for x in pal_factors(long_f[:3000]) if len(x) == L} for L in range(1, 31))
    check("B11c s_9 contains ALL palindromic factors of f of length <=30 (what the table needs)", ok)
    # table of letter counts
    expected = {0: (1, 0), 1: (1, 1), 2: (2, 1), 3: (3, 2), 4: (5, 3), 5: (8, 5), 6: (13, 8), 7: (21, 13),
                8: (34, 21), 9: (55, 34), 13: (377, 233), 14: (610, 377), 18: (4181, 2584), 19: (6765, 4181)}
    ok = all((s_word(k).count("0"), s_word(k).count("1")) == v for k, v in expected.items())
    check("B12 table of N_0, N_1 in s_k", ok)
    d20 = Decimal(F[20]) / Decimal(F[22])
    err = d20 - 1 / PHI ** 2
    pred = -(PSI ** 20) / (PHI ** 2 * F[22])
    check("B13 s_20: lengths 17711/10946/6765", (F[22], F[21], F[20]) == (17711, 10946, 6765))
    check("B14 density of 1 in s_20 = 0.3819660098 (thesis) and error -1.4e-9",
          f"{d20:.10f}" == "0.3819660098" and abs(float(err) + 1.4e-9) < 1e-10 and abs(err - pred) < Decimal("1e-40"),
          f"{d20:.10f}, err={float(err):.3e}")
    # error formula for all n
    ok = all(Decimal(F[n]) / Decimal(F[n + 2]) - 1 / PHI ** 2 + PSI ** n / (PHI ** 2 * F[n + 2]) < Decimal("1e-60")
             and abs(Decimal(F[n]) / Decimal(F[n + 2]) - 1 / PHI ** 2 + PSI ** n / (PHI ** 2 * F[n + 2])) < Decimal("1e-60")
             for n in range(1, 120))
    check("B15 exact error formula D(1,s_n)-1/phi^2 = -psi^n/(phi^2 F_{n+2})", ok)
    # incidence matrix
    M = np.array([[1, 1], [1, 0]], dtype=object)
    v = np.array([1, 0], dtype=object)
    ok = True
    x = v.copy()
    for n in range(1, 25):
        x = M.dot(x)
        ok &= (x[0] == F[n + 1] and x[1] == F[n])
    check("B16 (|s_n|_0,|s_n|_1)^T = M^n (1,0)^T", ok)


# ----------------------------------------------------------------------------
# C  Beatty
# ----------------------------------------------------------------------------
def sec_C():
    print("\n== C. Beatty description (Lemma mechanical) ==")
    f = s_word(23)
    N = 20000
    B1 = {fl_phi(m) for m in range(1, N)}
    B2 = {fl_phi2(m) for m in range(1, N)}
    ones = {n for n in range(1, N) if f[n - 1] == "1"}
    zeros = {n for n in range(1, N) if f[n - 1] == "0"}
    check("C1 positions of 1 = {floor(m phi^2)}", ones == {x for x in B2 if x < N})
    check("C2 positions of 0 = {floor(m phi)}", zeros == {x for x in B1 if x < N})
    check("C3 B1, B2 partition N>=1", not (B1 & B2) and set(range(1, 10000)) <= (B1 | B2))
    check("C4 B1 starts 1,3,4,6,8,9 ; B2 starts 2,5,7,10,13",
          [fl_phi(m) for m in range(1, 7)] == [1, 3, 4, 6, 8, 9] and [fl_phi2(m) for m in range(1, 6)] == [2, 5, 7, 10, 13])
    g1 = [fl_phi(n) for n in range(1, 31)]
    g2 = [fl_phi2(n) for n in range(1, 31)]
    check("C5 plotted Beatty data (g1(30)=48, g2(30)=78)", g1[-1] == 48 and g2[-1] == 78)


# ----------------------------------------------------------------------------
# D  framed words
# ----------------------------------------------------------------------------
def u_word(n):
    a, b = "b", "a"
    for _ in range(n):
        a, b = b, b + a
    return b


def sec_D():
    print("\n== D. Framed words q_m (Theorem thm1Asymptotic) ==")
    ok = u_word(0) == "a" and u_word(1) == "ab" and u_word(2) == "aba" and u_word(3) == "abaab"
    check("D1 u_0..u_3", ok)
    check("D2 q_1=aabb, q_2=aabab, q_3=aabaabb",
          ["a" + u_word(m) + "b" for m in (1, 2, 3)] == ["aabb", "aabab", "aabaabb"])
    ok1 = ok2 = ok3 = True
    for m in range(1, 25):
        q = "a" + u_word(m) + "b"
        ok1 &= (len(q) == F[m + 2] + 2 and q.count("a") == F[m + 1] + 1 and q.count("b") == F[m] + 1)
        d = Fraction(q.count("a"), len(q))
        ok2 &= (d == Fraction(F[m + 1] + 1, F[m + 2] + 2))
        dec = Decimal(d.numerator) / Decimal(d.denominator)
        formula = 1 / PHI - (1 / PHI ** 2 + PSI ** (m + 1)) / (PHI * (F[m + 2] + 2))
        ok3 &= abs(dec - formula) < Decimal("1e-60") and dec < 1 / PHI
    check("D3 |q_m|, |q_m|_a, |q_m|_b, m<=24 (brute force on strings)", ok1)
    check("D4 dens_a(q_m) = (F_{m+1}+1)/(F_{m+2}+2)", ok2)
    check("D5 exact error formula and dens_a(q_m)<1/phi, m<=24", ok3)
    table = {3: (0.571429, 0.428571, 0.600000, 0.400000), 5: (0.600000, 0.400000, 0.615385, 0.384615),
             7: (0.611111, 0.388889, 0.617647, 0.382353), 9: (0.615385, 0.384615, 0.617978, 0.382022),
             11: (0.617021, 0.382979, 0.618026, 0.381974), 13: (0.617647, 0.382353, 0.618033, 0.381967)}
    ok = True
    for m, row in table.items():
        q = "a" + u_word(m) + "b"
        u = u_word(m)
        vals = (q.count("a") / len(q), q.count("b") / len(q), u.count("a") / len(u), u.count("b") / len(u))
        ok &= all(abs(a - b) < 5e-7 for a, b in zip(vals, row))
    check("D6 Table (framed words) reproduced to 6 decimals", ok)
    ok = all(abs((Decimal(F[m + 1]) / F[m + 2]) - 1 / PHI + PSI ** (m + 1) / (PHI * F[m + 2])) < Decimal("1e-60")
             for m in range(1, 80))
    check("D7 dens_a(u_m)-1/phi = -psi^{m+1}/(phi F_{m+2})", ok)
    # Lemma fibA
    ok = all(abs(F[n + 1] - PHI * F[n] - PSI ** n) < Decimal("1e-30") for n in range(0, 120))
    check("D8 F_{n+1} - phi F_n = psi^n", ok)


# ----------------------------------------------------------------------------
# E  Millin and Fibonacci identities
# ----------------------------------------------------------------------------
def sec_E():
    print("\n== E. Telescoping series and identities ==")
    ok = all(Fraction(F[n - 1], F[n]) - Fraction(F[2 * n - 1], F[2 * n]) == Fraction((-1) ** n, F[2 * n]) for n in range(1, 40))
    check("E1 Lemma telesc: F_{n-1}/F_n - F_{2n-1}/F_{2n} = (-1)^n/F_{2n}", ok)
    ok = True
    for m in range(1, 9):
        for K in range(0, 5):
            M = m * 2 ** (K + 1)
            lhs = sum(Fraction(1, F[m * 2 ** k]) for k in range(1, K + 2))
            rhs = Fraction(F[m - 1], F[m]) - Fraction(F[M - 1], F[M]) + Fraction(1 - (-1) ** m, F[2 * m])
            ok &= (lhs == rhs)
    check("E2 Theorem millin finite identity (m<=8, K<=4), exact rational", ok)
    # infinite sum (m=1 gives (7-sqrt5)/2 including the term k=0)
    s = sum(Decimal(1) / F[2 ** k] for k in range(0, 9))  # F_{2^8}=F_256 enough
    check("E3 Millin: sum_{k>=0} 1/F_{2^k} = (7-sqrt5)/2", abs(s - (7 - Decimal(5).sqrt()) / 2) < Decimal("1e-50"))
    for m in (1, 2, 3, 4, 5):
        s = sum(Decimal(1) / F[m * 2 ** k] for k in range(1, 8))
        r = Decimal(F[m - 1]) / F[m] - 1 / PHI + Decimal(1 - (-1) ** m) / F[2 * m]
        check(f"E4 infinite formula m={m}", abs(s - r) < Decimal("1e-40"))
    # d'Ocagne
    ok = all(F[k + m] * F[k + 1] - F[k + m + 1] * F[k] == (-1) ** k * F[m] for k in range(1, 40) for m in range(0, 40))
    check("E5 d'Ocagne identity", ok)
    ok = True
    for a in range(1, 9):
        for b in range(1, 9):
            L = sum(Fraction((-1) ** k * F[a], F[k] * F[k + a]) for k in range(1, b + 1))
            R = sum(Fraction((-1) ** k * F[b], F[k] * F[k + b]) for k in range(1, a + 1))
            ok &= (L == R)
    check("E6 Prop. symm (a,b<=8), exact", ok)
    ok = True
    for r in range(1, 8):
        lhs = sum(Decimal(1) / (F[k] * F[k + 2 * r]) for k in range(1, 300))
        rhs = sum(Decimal(1) / (F[2 * n - 1] * F[2 * n]) for n in range(1, r + 1)) / F[2 * r]
        ok &= abs(lhs - rhs) < Decimal("1e-40")
    check("E7 Prop. sum2r (r<=7), 300 terms, error<1e-40", ok)


# ----------------------------------------------------------------------------
# F  palindromic prefixes
# ----------------------------------------------------------------------------
def sturm_char(ds, L):
    """s_{-1}=1, s_0=0, s_j = s_{j-1}^{d_j} s_{j-2}."""
    a, b = "1", "0"
    for d in ds:
        a, b = b, b * d + a
        if len(b) >= L:
            break
    return b[:L]


def pal_prefix_formula(ds, upto):
    q_prev, q = 1, 1  # q_{-1}, q_0
    out = set()
    k = 0
    while True:
        d_next = ds[k]
        for m in range(1, d_next + 1):
            L = m * q + q_prev - 2
            if L > 0 and L <= upto:
                out.add(L)
        q_prev, q = q, d_next * q + q_prev
        k += 1
        if k >= len(ds) or q_prev > upto + 5:
            break
    return sorted(out)


def sec_F():
    print("\n== F. Palindromic prefixes and palindrome density (Sect. 3.6) ==")
    f = s_word(22)
    pl = pal_prefix_lengths(f, 4000)
    check("F1 palindromic prefixes of f have lengths exactly F_{k+2}-2 (<=4000)",
          pl == sorted(x for x in {F[k + 2] - 2 for k in range(2, 30)} if 0 < x <= 4000), str(pl[:10]))
    ok = all(s_word(n)[:-2] == s_word(n)[:-2][::-1] for n in range(1, 20))
    check("F2 Lemma palFib: s_n minus last two letters is a palindrome", ok)
    # ternary word y
    ys = ["01", "02"]
    for n in range(2, 20):
        ys.append(ys[-1] + ys[-2])
    check("F3 y_2=0201, y_3=020102, y_4=0201020201", (ys[2], ys[3], ys[4]) == ("0201", "020102", "0201020201"))
    ok = all(len(ys[n]) == 2 * F[n + 1] and ys[n].count("0") == F[n + 1]
             and ys[n].count("1") == (F[n - 1] if n >= 1 else 1) and ys[n].count("2") == F[n] for n in range(1, 18))
    check("F4 Prop. lemmanovelwordn1: counts of y_n (n>=1)", ok)
    exp = {2: (4, 2, 1, 1), 5: (16, 8, 3, 5), 9: (110, 55, 21, 34), 15: (1974, 987, 377, 610)}
    check("F5 Table of |y_n| etc. (n=2,5,9,15)",
          all((len(ys[n]), ys[n].count("0"), ys[n].count("1"), ys[n].count("2")) == v for n, v in exp.items()))
    # y_inf = iota_0(rho(f))
    rho = {"0": "2", "1": "1"}
    yinf = "".join("0" + rho[c] for c in s_word(20))
    check("F6 y_inf = iota_0(rho(f)) (prefix check)", ys[18] == yinf[:len(ys[18])])
    pl = pal_prefix_lengths(yinf[:9000], 8000)
    check("F7 palindromic prefixes of y: 2F_k-3, k>=3 (1,3,7,13,23,39,65,...)",
          pl == sorted(x for x in {2 * F[k] - 3 for k in range(3, 40)} if x <= 8000), str(pl[:8]))
    # richness
    ok = True
    for n in range(0, 15):
        ok &= (len(pal_factors(ys[n])) == len(ys[n]))
    check("F8 Theorem Yrich: y_n is rich (|P(y_n)|=|y_n|), n<=14", ok)
    z = ys[10] + "0"
    check("F9 Lemma insertcount: |P(z)|=1+2|P(x)|",
          len(pal_factors(z)) == 1 + 2 * len(pal_factors("".join(c for c in z if c != "0"))))
    # check |P(z)| = |z| for z = y_n 0
    check("F10 |P(y_n 0)| = |y_n 0| (n<=12)", all(len(pal_factors(ys[n] + "0")) == len(ys[n]) + 1 for n in range(0, 13)))
    # proposition mk
    ok = True
    for k in range(4, 17):
        L = 2 * F[k] - 3
        mk = yinf[:L]
        ms = mk[1:-1]
        ok &= (mk == mk[::-1] and ms == ms[::-1])
        ok &= (mk.count("0"), mk.count("1"), mk.count("2")) == (F[k] - 1, F[k - 2] - 1, F[k - 1] - 1)
        ok &= (ms.count("0"), ms.count("1"), ms.count("2")) == (F[k] - 3, F[k - 2] - 1, F[k - 1] - 1)
    check("F11 Prop. mk: letter counts of m_k and m_k*", ok)

    # characteristic Sturmian words: Theorem ct:palprefSturm
    random.seed(1)
    ok = True
    tests = [[1] * 40, [1] + [2] * 40, [1] + [3] * 30, [2, 1, 3, 1, 2, 2, 1, 4, 1, 1, 2, 3, 1, 1, 5, 1, 2] + [1] * 20]
    for _ in range(4):
        tests.append([random.randint(1, 4) for _ in range(40)])
    for ds in tests:
        w = sturm_char(ds, 9000)
        got = pal_prefix_lengths(w, 4000)
        want = pal_prefix_formula(ds, 4000)
        ok &= (got == want)
    check("F12 palindromic prefix lengths = {m q_k + q_{k-1} - 2} for 8 slopes", ok)

    # Corollary palSturm (5): limits
    def delta_from_formula(ds, upto=10 ** 60):
        L = pal_prefix_formula(ds, upto)
        r = [Fraction(L[i + 1], L[i]) for i in range(len(L) - 1)]
        tail = r[len(r) // 2:]
        return float(max(tail))
    for a in (1, 2, 3, 4, 6):
        ds = [1] + [a] * 80
        d = delta_from_formula(ds)
        ta = (math.sqrt(a * a + 4) - a) / 2
        pred = 1 + 1 / (1 + ta)
        check(f"F13 delta for a_k={a}: limsup n_(i+1)/n_i ~ 1+1/(1+t_a)={pred:.5f}",
              abs(d - pred) < 2e-3, f"computed {d:.5f}; dens_p={1/pred:.5f}; (1+t)/(2+t)={(1+ta)/(2+ta):.5f}")
    check("F14 dens_p(a=2) = 2-sqrt2", abs((1 + (math.sqrt(8) - 2) / 2) / (2 + (math.sqrt(8) - 2) / 2) - (2 - math.sqrt(2))) < 1e-12)
    # unbounded partial quotients: ratio -> 2
    ds = [1] + [n for n in range(2, 40)] * 1
    d = delta_from_formula(ds)
    check("F15 growing partial quotients: ratio approaches 2 (dens_p -> 1/2)", d > 1.9, f"{d:.4f}")
    # Prop. dichotomy lower bound
    ok = True
    for A in range(1, 8):
        for _ in range(30):
            ds = [random.randint(1, A) for _ in range(60)]
            ds[0] = max(ds[0], 1)
            d = delta_from_formula(ds)
            ok &= (d <= (2 * A + 3) / (A + 2) + 1e-9)
    check("F16 Prop. dichotomy (2): delta <= (2A+3)/(A+2) on 210 random slopes", ok)
    # iota_0 invariance
    w = sturm_char([1, 2, 1, 3, 1, 1, 2] + [1] * 30, 3000)
    x = "".join("0" + {"0": "2", "1": "1"}[c] for c in w[:1500])
    L1 = pal_prefix_lengths(w, 1400)
    L2 = pal_prefix_lengths(x, 2800)
    check("F17 Lemma interleave: pal. prefixes of iota_0(x) are 1 and 2m+1", L2 == sorted({1} | {2 * m + 1 for m in L1 if 2 * m + 1 <= 2800}))

    # m-bonacci
    def mbonacci(m, L):
        mor = {chr(97 + i): "a" + chr(98 + i) for i in range(m - 1)}
        mor[chr(97 + m - 1)] = "a"
        w = "a"
        while len(w) < L:
            w = "".join(mor[c] for c in w)
        return w
    for m, beta in ((3, 1.839286755214161), (4, 1.9275619754829254)):
        w = mbonacci(m, 40000)
        L = pal_prefix_lengths(w, 15000)
        r = [L[i + 1] / L[i] for i in range(len(L) - 1)]
        tail = r[len(r) // 2:]
        # recurrence n_{i+1} = 2 n_i - n_{i-m}
        ok_rec = False
        for off in (0, 1):
            LL = ([0] if off else []) + L
            ok_rec |= all(LL[i + 1] == 2 * LL[i] - LL[i - m] for i in range(m + 1, len(LL) - 1))
        check(f"F18 m={m}: palindromic-prefix ratio tail max ~ beta_m={beta:.4f}", abs(max(tail) - beta) < 0.08,
              f"max tail ratio={max(tail):.4f}, #pal.prefixes={len(L)}")
        check(f"F19 m={m}: recurrence n_(i+1)=2n_i-n_(i-m) holds (i>m)", ok_rec)
    # beta_m is the root
    ok = True
    for m, beta in ((3, 1.839286755214161), (4, 1.9275619754829254)):
        ok &= abs(beta ** m - sum(beta ** i for i in range(m))) < 1e-9
    check("F20 beta_3=1.8392.., beta_4=1.9276.. are roots of x^m = x^{m-1}+...+1", ok)


# ----------------------------------------------------------------------------
# G  Generalized Thue-Morse
# ----------------------------------------------------------------------------
def digit_sum_array(L, b):
    s = np.zeros(L, dtype=np.int64)
    x = np.arange(L, dtype=np.int64)
    while x.any():
        s += x % b
        x //= b
    return s


def longest_run(mask):
    m = np.concatenate(([0], mask.astype(np.int8), [0]))
    d = np.diff(m)
    st = np.where(d == 1)[0]
    en = np.where(d == -1)[0]
    return int((en - st).max()) if len(st) else 0


def sec_G():
    print("\n== G. Generalized Thue-Morse words W^(b) (Sect. 3.8) ==")
    # Prop gTM: fixed point of xi_b equals s_b(n) mod N
    ok = True
    for N in range(2, 7):
        for b in range(2, 7):
            w = [0]
            while len(w) < 5000:
                w = [(t + r) % N for t in w for r in range(b)]
            s = digit_sum_array(len(w), b) % N
            ok &= bool((np.array(w) == s).all())
    check("G1 Prop. gTM: fixed point of xi_b = s_b(n) mod N (N,b<=6)", ok)
    # example N=b=3
    s = digit_sum_array(100, 3) % 3
    check("G2 Example N=b=3: W[24..29] = 1,2,0,1,2,0", list(s[24:30]) == [1, 2, 0, 1, 2, 0])
    # kappa
    bad = []
    n_tested = 0
    for N in range(2, 7):
        for b in range(2, 7):
            d = math.gcd(N, b - 1)
            g = N // d
            if g == 1:
                continue
            L = min(300000, max(30000, 4 * b ** g))
            s = digit_sum_array(L + 200, b)
            for j in range(N):
                jp = j if j != 0 else N
                kappa = max(2, (2 * b) // jp)
                best = 0
                for m in range(1, 101):
                    c = s[m:L + m] - s[:L]
                    run = longest_run((c % N) == (j % N))
                    if run >= m:
                        best = max(best, run // m + 1)
                n_tested += 1
                if best != kappa:
                    bad.append((N, b, j, best, kappa))
    check(f"G3 Theorem tw:main: kappa(N,b,j) equals brute-force max k, {n_tested} triples (m<=100)",
          not bad, f"mismatches (N,b,j,found,formula): {bad[:6]}")
    # b=2 corollary
    ok = True
    for N in range(2, 9):
        s = digit_sum_array(120000, 2)
        for j in range(N):
            best = 0
            for m in range(1, 200):
                c = s[m:100000 + m] - s[:100000]
                run = longest_run((c % N) == j)
                if run >= m:
                    best = max(best, run // m + 1)
            ok &= (best == (4 if j == 1 % N else 2))
    check("G4 Corollary: b=2: kappa=4 if j=1 mod N else 2 (N<=8, m<200)", ok)
    # N=2,b=2 : strong (4,sigma) repetitions abab exist, none with k=5
    t = "".join(map(str, digit_sum_array(2 ** 16, 2) % 2))
    check("G5 Thue-Morse contains 0101 and 1010 (strong (4,sigma)) but not 01010 / 10101 (k=5)",
          "0101" in t and "1010" in t and "01010" not in t and "10101" not in t)
    # Theorem avoidance witness i=2^N-2
    ok = True
    for N in range(2, 9):
        i = 2 ** N - 2
        s = digit_sum_array(i + 5, 2) % N
        ok &= (s[i + 1] - s[i]) % N == 1 and (s[i + 2] - s[i + 1]) % N == 1
    check("G6 witness W[i..i+2] at i=2^N-2 is a strong (3,sigma)-repetition (N<=8)", ok)
    # Observation: complexity
    def complexity(N, size, K=40):
        w = bytes((digit_sum_array(size, 2) % N).astype(np.uint8))
        return [len({w[i:i + k] for i in range(len(w) - k + 1)}) for k in range(1, K + 1)]
    for N in (2, 3, 4):
        p1 = complexity(N, 2 ** 16)
        p2 = complexity(N, 2 ** 17)
        diffs = {p2[k] - p2[k - 1] for k in range(1, 40)}
        check(f"G7 N={N}: p(k) stable (2^16 vs 2^17), p(1)={N}, p(2)={N*N}, differences in {{N(N-1),N^2}}",
              p1 == p2 and p2[0] == N and p2[1] == N * N and diffs <= {N * (N - 1), N * N},
              f"p(40)={p2[39]}, differences={sorted(diffs)}")
        if N == 3:
            check("G8 N=3: p(40)=303 (thesis)", p2[39] == 303, f"computed {p2[39]}")
    # Dynamics: letter-frequency deviation (Theorem dynWb (4)): sanity
    b, N = 3, 2
    s = digit_sum_array(3 ** 12, b) % N
    dev = max(abs(int((s[i:i + 2000] == 1).sum()) - 1000) for i in range(0, 3 ** 12 - 2000, 997))
    check("G9 Thm dynWb(4): window deviation from n/N is o(n) (b=3,N=2, n=2000)", dev < 100, f"max deviation {dev}")


# ----------------------------------------------------------------------------
# H  factor trees
# ----------------------------------------------------------------------------
def tree_stats(w, n):
    V = {""}
    for k in range(1, n + 1):
        V |= {w[i:i + k] for i in range(len(w) - k + 1)}
    deg = {}
    for v in V:
        ch = sum(1 for c in "01" if (v + c) in V) if len(v) < n else 0
        deg[v] = ch + (1 if v else 0)
    edges = [(v[:-1], v) for v in V if v]
    irr = sum(abs(deg[a] - deg[b]) for a, b in edges)
    irr2 = sum(abs(deg[a] ** 2 - deg[b] ** 2) for a, b in edges)
    M1 = sum(x * x for x in deg.values())
    M2 = sum(deg[a] * deg[b] for a, b in edges)
    return len(V), irr, irr2, M1, M2


def sec_H():
    print("\n== H. Factor trees of Sturmian words (Theorem factortree) ==")
    random.seed(7)
    ok = True
    cnt = 0
    fails = []
    for d in range(1, 6):
        for trial in range(4):
            ds = [d] + [random.randint(1, 3) for _ in range(40)]
            if trial == 0:
                ds = [d] + [1] * 40
            w = sturm_char(ds, 30000)
            assert w.startswith("0" * d + "1")
            for n in range(d + 2, 15):
                V, irr, irr2, M1, M2 = tree_stats(w, n)
                exp = ((n + 1) * (n + 2) // 2, 4 * n - 2 * d, 18 * n - 10 * d - 2, 2 * n * n + 8 * n - 4,
                       2 * n * n + 10 * n - 11 + d)
                cnt += 1
                if (V, irr, irr2, M1, M2) != exp:
                    ok = False
                    fails.append((d, trial, n, (V, irr, irr2, M1, M2), exp))
    check(f"H1 |V|, irr, irr_2, M_1, M_2 formulas ({cnt} (word,n) pairs, d=1..5)", ok, str(fails[:2]))
    rows = []
    ok = True
    expected_table = {1: (45, 30, 132, 188, 198), 2: (45, 28, 122, 188, 199), 3: (45, 26, 112, 188, 200), 4: (45, 24, 102, 188, 201)}
    for d in range(1, 5):
        w = sturm_char([d] + [1] * 40, 30000)
        st = tree_stats(w, 8)
        ok &= (st == expected_table[d])
        rows.append((d, st))
    check("H2 Table T_8 for d=1..4 reproduced", ok, str(rows[0]))
    # Fibonacci: irr=4n-2, irr2=18n-12, M2=2n^2+10n-10
    f = s_word(22)
    ok = all((lambda st, n: st[1] == 4 * n - 2 and st[2] == 18 * n - 12 and st[4] == 2 * n * n + 10 * n - 10)(tree_stats(f, n), n)
             for n in range(3, 16))
    check("H3 Fibonacci word: irr=4n-2, irr_2=18n-12, M_2=2n^2+10n-10", ok)
    # P4 values
    P4 = nx.path_graph(4)
    A = sum(abs(P4.degree(u) - P4.degree(v)) for u, v in P4.edges())
    A1 = sum(Fraction((P4.degree(u) - P4.degree(v)) ** 2, nx.shortest_path_length(P4, u, v)) for u, v in
             [(a, b) for a in P4 for b in P4 if a < b])
    G = sum(P4.degree(u) * P4.degree(v) * nx.shortest_path_length(P4, u, v) for u, v in [(a, b) for a in P4 for b in P4 if a < b])
    check("H4 A(P4)=2, A_1(P4)=3, Gut(P4)=19", (A, A1, G) == (2, 3, 19), f"{A},{A1},{G}")


# ----------------------------------------------------------------------------
# I  trees
# ----------------------------------------------------------------------------
def irr_idx(G):
    return sum(abs(G.degree(u) - G.degree(v)) for u, v in G.edges())


def irr2_idx(G):
    return sum(abs(G.degree(u) ** 2 - G.degree(v) ** 2) for u, v in G.edges())


def sec_I():
    print("\n== I. Trees: Prop. pron1, conjecture, Lemma lemn01, Albertson bounds ==")
    ok1 = ok2 = ok3 = ok4 = True
    cnt = 0
    bad = []
    for n in range(5, 15):
        for T in nx.nonisomorphic_trees(n):
            cnt += 1
            D = max(dict(T.degree()).values())
            r2 = irr2_idx(T)
            big = sum(1 for v in T if T.degree(v) > 2)
            path = D <= 2
            ok1 &= r2 >= D * (D * D - 1)
            eq = r2 == D * (D * D - 1)
            ok2 &= (eq == (path or big == 1))
            star = D == n - 1
            if not star:
                if not (r2 < n * (n - 1) * (n - 2)):
                    ok3 = False
                    bad.append(n)
            A = irr_idx(T)
            leaves = sum(1 for v in T if T.degree(v) == 1)
            ok4 &= (A >= leaves and A % 2 == 0 and ((A == 2) == path))
    check(f"I1 Prop. pron1: irr_2 >= D(D^2-1) on all {cnt} trees with 5<=n<=14", ok1)
    check("I2 equality iff path or exactly one vertex of degree>2", ok2)
    check("I3 Conjecture irr_2(T) < n(n-1)(n-2) for all non-star trees, n=5..14", ok3, f"counterexamples n={bad[:3]}")
    check("I4 Prop. Albertson: A(T)>=leaves, A even, A=2 iff path (n=5..14)", ok4)
    # Lemma lemn01 on random graphs
    random.seed(3)
    ok = True
    trials = 0
    for _ in range(3000):
        n = random.randint(4, 10)
        G = nx.gnp_random_graph(n, random.uniform(0.2, 0.7), seed=random.randint(0, 10 ** 9))
        non = [(u, v) for u in G for v in G if u < v and not G.has_edge(u, v)]
        if not non:
            continue
        a, b = random.choice(non)
        if G.degree(a) < G.degree(b):
            a, b = b, a
        da, db = G.degree(a), G.degree(b)
        ra = sum(1 for x in G[a] if G.degree(x) > da)
        rb = sum(1 for x in G[b] if G.degree(x) > db)
        before = irr2_idx(G)
        G.add_edge(a, b)
        after = irr2_idx(G)
        pred = 3 * da * (da + 1) + db * (db - 1) - 2 * ((2 * da + 1) * ra + (2 * db + 1) * rb)
        trials += 1
        ok &= (after - before == pred)
    check(f"I5 Lemma lemn01 verified on {trials} random graph/non-edge pairs", ok)
    # A_0 identity, sigma limit
    ok = True
    for n in range(4, 11):
        for T in nx.nonisomorphic_trees(n):
            deg = dict(T.degree())
            A0 = sum((deg[u] - deg[v]) ** 2 for u in T for v in T if u < v)
            M1 = sum(x * x for x in deg.values())
            ok &= (A0 == n * M1 - 4 * (n - 1) ** 2)
    check("I6 A_0(G) = n M_1 - 4 m^2 (all trees n<=10)", ok)
    # Gut_0
    ok = True
    for n in range(4, 10):
        for T in nx.nonisomorphic_trees(n):
            deg = dict(T.degree())
            g0 = sum(deg[u] * deg[v] for u in T for v in T if u < v)
            M1 = sum(x * x for x in deg.values())
            ok &= (2 * g0 == 4 * (n - 1) ** 2 - M1)
    check("I7 Gut_0(G) = (4m^2 - M_1)/2 (all trees n<=9)", ok)
    star = nx.star_graph(7)
    check("I8 star S_n: irr_2 = n(n-1)(n-2), irr=(n-1)(n-2)", irr2_idx(star) == 8 * 7 * 6 and irr_idx(star) == 7 * 6)


# ----------------------------------------------------------------------------
# J  2D recursion
# ----------------------------------------------------------------------------
def sec_J():
    print("\n== J. Two-dimensional recursion T(u,v) ==")
    memo = {}

    def T(u, v):
        if (u, v) in memo:
            return memo[(u, v)]
        if u == 0:
            r = "a"
        elif v == 0:
            r = "b"
        else:
            r = T(u - 1, v) + T(u, v - 1)
        memo[(u, v)] = r
        return r
    check("J1 T(1,1)=ab, T(2,1)=abb, T(1,2)=aab, T(2,2)=aababb",
          (T(1, 1), T(2, 1), T(1, 2), T(2, 2)) == ("ab", "abb", "aab", "aababb"))
    ok = True
    for u in range(1, 11):
        for v in range(1, 11):
            w = T(u, v)
            ok &= (len(w) == math.comb(u + v, u) and w.count("a") == math.comb(u + v - 1, u)
                   and w.count("b") == math.comb(u + v - 1, v)
                   and Fraction(w.count("a"), len(w)) == Fraction(v, u + v))
    check("J2 Prop. newprove: lengths, letter counts, densities v/(u+v), u,v<=10", ok)
    check("J3 lattice grid in Fig.: central cells 2,6,20,70", [math.comb(2 * i, i) for i in range(1, 5)] == [2, 6, 20, 70])


# ----------------------------------------------------------------------------
# K  factorial word
# ----------------------------------------------------------------------------
def sec_K():
    print("\n== K. Word obtained by concatenating base-b expansions of n! (finite sanity check) ==")
    for b, fmt, Nmax, kmax in ((2, lambda x: bin(x)[2:], 1500, 7), (10, str, 3000, 3), (16, lambda x: format(x, "x"), 1500, 2)):
        parts = []
        f = 1
        for n in range(0, Nmax + 1):
            if n > 0:
                f *= n
            parts.append(fmt(f))
        S = "".join(parts)
        digits = "0123456789abcdef"[:b]
        covered = 0
        import itertools
        for k in range(1, kmax + 1):
            if all("".join(t) in S for t in itertools.product(digits, repeat=k)):
                covered = k
            else:
                break
        check(f"K1 base {b}: all words of length k<={kmax} occur in the concatenation of 0!..{Nmax}!",
              covered == kmax, f"covered up to k={covered}, prefix length {len(S)}")
    check("K2 F_10 starts 112624120720504040320", "".join(str(math.factorial(n)) for n in range(9)) == "112624120720504040320")


# ----------------------------------------------------------------------------
# L  statistics of chapter 4 from reported summary numbers
# ----------------------------------------------------------------------------
def sec_L():
    print("\n== L. Chapter 4: recomputation from the summary statistics printed in the thesis ==")
    from scipy import stats
    n = 42
    T = {  # (mean, sd) pre A, pre B, post A, post B
        "think": ((2.71, 2.30), (3.76, 2.18), (8.26, 2.87), (9.26, 2.65), (5.55, 3.62), (5.50, 3.58)),
        "nn": ((47.79, 24.68), (59.74, 22.11), (65.36, 21.61), (77.83, 20.19), (17.57, 21.08), (18.10, 19.50)),
    }

    def welch(a, b):
        r = stats.ttest_ind_from_stats(b[0], b[1], n, a[0], a[1], n, equal_var=False)
        se = math.sqrt(a[1] ** 2 / n + b[1] ** 2 / n)
        df = (a[1] ** 2 / n + b[1] ** 2 / n) ** 2 / ((a[1] ** 2 / n) ** 2 / (n - 1) + (b[1] ** 2 / n) ** 2 / (n - 1))
        tcrit = stats.t.ppf(0.975, df)
        diff = b[0] - a[0]
        sp = math.sqrt((a[1] ** 2 + b[1] ** 2) / 2)
        return diff, (diff - tcrit * se, diff + tcrit * se), r.statistic, df, r.pvalue, diff / sp

    exp = {("think", "pre"): (1.05, 2.14, 0.035, 0.47), ("think", "post"): (1.00, 1.66, 0.101, 0.36),
           ("nn", "pre"): (11.95, 2.34, 0.022, 0.51), ("nn", "post"): (12.48, 2.73, 0.0077, 0.60),
           ("think", "gain"): (-0.05, -0.06, 0.95, None), ("nn", "gain"): (0.52, 0.12, 0.91, None)}
    idx = {"pre": (0, 1), "post": (2, 3), "gain": (4, 5)}
    for (test, stage), (ed, et, ep, ecd) in exp.items():
        a, b = T[test][idx[stage][0]], T[test][idx[stage][1]]
        diff, ci, t, df, p, d = welch(a, b)
        good = abs(diff - ed) < 0.011 and abs(t - et) < 0.02 and abs(p - ep) < max(0.002, 0.05 * ep)
        if ecd is not None:
            good &= abs(d - ecd) < 0.011
        check(f"L1 Welch {test}/{stage}: diff={diff:.2f}, t={t:.2f}, p={p:.4f}, CI=({ci[0]:.2f},{ci[1]:.2f}), d={d:.2f}", good)
    # reliability CI for r
    r, N = 0.855, 40
    z = math.atanh(r)
    se = 1 / math.sqrt(N - 3)
    lo, hi = math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)
    check("L2 pilot reliability: 95% CI for r=0.855, N=40 is (0.740, 0.921)", abs(lo - 0.740) < 0.002 and abs(hi - 0.921) < 0.002,
          f"({lo:.3f},{hi:.3f})")
    # paired shift: 9.18 vs 8.15, t(39)=4.29 -> sd_diff
    sd_diff = 1.03 * math.sqrt(40) / 4.29
    check("L3 pilot paired shift t(39)=4.29 implies p<0.001", 2 * stats.t.sf(4.29, 39) < 0.001, f"sd of differences ~{sd_diff:.2f}")
    # sensitivity table (noncentral t)
    def detectable(sd, power, n=42, alpha=0.05):
        df = 2 * n - 3
        tc = stats.t.ppf(1 - alpha / 2, df)
        lo_, hi_ = 0.0, 10.0 * sd
        for _ in range(100):
            mid = (lo_ + hi_) / 2
            ncp = mid / (sd * math.sqrt(2 / n))
            pw = stats.nct.sf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp)
            if pw < power:
                lo_ = mid
            else:
                hi_ = mid
        return hi_
    vals = {(0.8, 2.78): 1.72, (0.9, 2.78): 1.99, (0.8, 17.05): 10.55, (0.9, 17.05): 12.21}
    ok = True
    out = []
    for (pw, sd), e in vals.items():
        v = detectable(sd, pw)
        out.append(f"{v:.2f}")
        ok &= abs(v - e) < 0.03
    check("L4 Table tab4 detectable differences (80%/90% power; sd 2.78 and 17.05)", ok, "computed " + ", ".join(out))
    # adjusted effect size claims: adjusted d ~ 0.37
    check("L5 adjusted d ~ 1.03/2.78 = 0.37 and 6.23/17.05 = 0.37", abs(1.03 / 2.78 - 0.37) < 0.005 and abs(6.23 / 17.05 - 0.37) < 0.005)
    # partial eta^2 from t: t^2/(t^2+df)
    e1 = 1.65 ** 2 / (1.65 ** 2 + 81)
    e2 = 1.62 ** 2 / (1.62 ** 2 + 81)
    check("L6 partial eta^2 from t(81): 0.033 and 0.031", abs(e1 - 0.033) < 0.001 and abs(e2 - 0.031) < 0.001, f"{e1:.4f}, {e2:.4f}")
    p1 = 2 * stats.t.sf(1.65, 81)
    p2 = 2 * stats.t.sf(1.62, 81)
    check("L7 ANCOVA p-values from t(81): 0.102 and 0.109", abs(p1 - 0.102) < 0.002 and abs(p2 - 0.109) < 0.002, f"{p1:.3f}, {p2:.3f}")
    print("   NOTE: ANCOVA, paired t-tests and Mann-Whitney need the raw data; see README (data/ folder).")


SECTIONS = {"A": sec_A, "B": sec_B, "C": sec_C, "D": sec_D, "E": sec_E, "F": sec_F, "G": sec_G,
            "H": sec_H, "I": sec_I, "J": sec_J, "K": sec_K, "L": sec_L}

if __name__ == "__main__":
    which = [a.upper() for a in sys.argv[1:]] or list(SECTIONS)
    for k in which:
        SECTIONS[k]()
    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n==== {len(RESULTS) - len(bad)}/{len(RESULTS)} checks passed ====")
    for n in bad:
        print("FAILED:", n)
    sys.exit(1 if bad else 0)
