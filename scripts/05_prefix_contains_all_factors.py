"""Which prefixes s_n of the Fibonacci word contain ALL factors of length <=30?  Palindromes per length."""
def fib_word(n):
    a, b = "1", "0"
    for _ in range(n):
        a, b = b, b + a
    return b

def factors(w, L):
    return {w[i:i + L] for i in range(len(w) - L + 1)}

def pals(S):
    return {u for u in S if u == u[::-1]}

if __name__ == "__main__":
    big = fib_word(25)                       # length 196418
    full = {L: factors(big, L) for L in range(1, 31)}
    print("n, |s_n|, contains all factors of length<=30?")
    for n in range(5, 16):
        s = fib_word(n)
        print(n, len(s), all(factors(s, L) == full[L] for L in range(1, 31)))
    print("L, number of palindromic factors (expect 2 for odd L, 1 for even L)")
    for L in range(1, 31):
        print(L, len(pals(full[L])))
