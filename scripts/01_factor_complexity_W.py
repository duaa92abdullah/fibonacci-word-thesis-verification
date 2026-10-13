"""Factor complexity of the generalized Thue-Morse word W (Observation in Sect. 3.8).
Expected: p(1)=N, p(2)=N^2, differences p(k+1)-p(k) in {N(N-1), N^2}; N=3: p(40)=303."""
def W(N, b, L):
    def sb(n):
        s = 0
        while n:
            s += n % b
            n //= b
        return s
    return [sb(n) % N for n in range(L)]

def p(w, k):
    return len({tuple(w[i:i + k]) for i in range(len(w) - k + 1)})

if __name__ == "__main__":
    for N in (2, 3, 4):
        w = W(N, 2, 100000)
        pk = [p(w, k) for k in range(1, 41)]
        diffs = sorted({pk[i + 1] - pk[i] for i in range(len(pk) - 1)})
        print(f"N={N}: p(1..40)={pk}\n     differences={diffs}  (expected subset of {{{N*(N-1)}, {N*N}}})")
