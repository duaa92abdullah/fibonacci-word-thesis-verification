"""Figure for the thesis (figfib_conv.pdf): density of the letter 1 in s_n -> 1/phi^2 (damped oscillation).
Uses the thesis convention s_n = Theta^n(0), |s_n|_1 = F_n, |s_n| = F_{n+2}."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

F = [0, 1]
for _ in range(40):
    F.append(F[-1] + F[-2])
phi = (1 + 5 ** 0.5) / 2
ns = list(range(1, 23))
rho1 = [F[n] / F[n + 2] for n in ns]        # density of 1 in s_n

plt.figure(figsize=(9, 5))
plt.plot(ns, rho1, "o-", color="#1f77b4", markersize=5, label=r"Empirical density of $1$ in $s_n$")
plt.axhline(1 / phi ** 2, color="#d62728", linestyle="--", linewidth=1.5, label=r"$\varphi^{-2}\approx0.3819$")
plt.title(r"Convergence of the frequency of $1$ in $s_n$"); plt.xlabel(r"Index $n$"); plt.ylabel(r"$|s_n|_1/|s_n|$")
plt.xticks(range(1, 23, 2)); plt.grid(True, linestyle=":", alpha=0.6); plt.legend(loc="upper right")
plt.savefig("figfib_conv.pdf", bbox_inches="tight")
print([round(x, 6) for x in rho1[:8]], "... saved figfib_conv.pdf")
