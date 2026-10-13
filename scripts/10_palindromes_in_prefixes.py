"""Density of palindromic factors in prefixes s_n (thesis convention s_n = Theta^n(0), so the rare letter is 1).
Reproduces the table:  00 in s_4: 1/8 ; 101 in s_5: 1/13 ; 010 in s_5: 4/13 ; 1001 in s_6: 4/21.
(The original script used the complemented words f_n, which swap 0<->1 and therefore print 11, 010, 101, 0110.)"""
def fib_word(n):
    a, b = "1", "0"
    for _ in range(n):
        a, b = b, b + a
    return b

def analyze_palindromes(word, length):
    counts = {}
    for i in range(len(word) - length + 1):
        sub = word[i:i + length]
        if sub == sub[::-1]:
            counts[sub] = counts.get(sub, 0) + 1
    return {k: f"{v}/{len(word)}" for k, v in counts.items()}

if __name__ == "__main__":
    for n, L in ((4, 2), (5, 3), (6, 4)):
        w = fib_word(n)
        print(f"s_{n} = {w} (length {len(w)}): palindromes of length {L}: {analyze_palindromes(w, L)}")
