"""Convergence of the density of '101' in prefixes of the Fibonacci word (limit 1/phi^4)."""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def fib_prefix(n):
    a, b = "0", "01"
    while len(b) < n:
        a, b = b, b + a
    return b[:n]

def count_occurrences(word, pattern):
    m = len(pattern)
    return sum(1 for i in range(len(word) - m + 1) if word[i:i + m] == pattern)

if __name__ == "__main__":
    phi = (1 + math.sqrt(5)) / 2
    lengths = [200, 400, 600, 800, 1000, 1500, 2000]
    dens = [count_occurrences(fib_prefix(n), "101") / n for n in lengths]
    plt.figure()
    plt.plot(lengths, dens, marker="o", label="density of 101")
    plt.axhline(1 / phi ** 4, color="r", linestyle="--", label=r"$1/\varphi^4\approx0.1459$")
    plt.xlabel("Prefix length N"); plt.ylabel("Density of 101"); plt.legend()
    plt.title("Convergence of density in the Fibonacci word")
    plt.savefig("density_101.pdf", bbox_inches="tight")
    print(list(zip(lengths, [round(d, 6) for d in dens])), "-> saved density_101.pdf")
