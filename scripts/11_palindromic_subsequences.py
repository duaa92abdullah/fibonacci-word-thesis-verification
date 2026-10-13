"""Palindromic factors versus scattered palindromic subsequences (Sect. 2.1.2)."""
def palindromic_substrings(s):
    n, out = len(s), set()
    for c in range(n):
        for l, r in ((c, c), (c, c + 1)):          # odd and even centres
            while l >= 0 and r < n and s[l] == s[r]:
                out.add(s[l:r + 1]); l -= 1; r += 1
    return out

def palindromic_subsequences(s):
    out = {""}
    for ch in s:                                   # all subsequences (2^n worst case), keep palindromes at the end
        out |= {x + ch for x in out}
    return {x for x in out if x and x == x[::-1]}

if __name__ == "__main__":
    for w in ("abaa", "abab", "abacaba"):
        P, SP = palindromic_substrings(w), palindromic_subsequences(w)
        print(f"{w}: P={len(P)} {sorted(P)}  SP={len(SP)} {sorted(SP)}")
    assert len(palindromic_substrings("abaa")) == 4 and len(palindromic_subsequences("abaa")) == 5
    assert len(palindromic_substrings("abab")) == 4 and len(palindromic_subsequences("abab")) == 6
    print("matches the examples in the thesis (abaa: 4 and 5; abab: 4 and 6)")
