from sage.all import *
from graph_plotting import compute_and_plot_gso, convert_to_fpylll
import random

def lattice_attack(ciphertexts, moduli, e):
    assert len(ciphertexts) == len(moduli), "The number of ciphertexts and N must be the same"

    random_pad_factors = [c[0] for c in ciphertexts]
    random_pad_offsets = [c[1] for c in ciphertexts]
    ciphertext = [c[2] for c in ciphertexts]
    ni = moduli  # use moduli list

    M = Matrix.identity(len(ciphertexts))
    ti = []
    for i in range(len(ciphertexts)):
        ti.append(crt(list(M[i]), ni))  # use the Chinese Remainder Theorem on M

    Zmodn = Zmod(prod(ni))
    R.<x> = PolynomialRing(Zmodn)

    gi = []
    for i in range(len(ciphertexts)):
        gi.append(ti[i] * ((random_pad_factors[i] * x + random_pad_offsets[i]) ** e - ciphertext[i]))

    g = sum(gi)
    g = g.monic()

    roots = adapted_small_roots(g)

    if roots:
        return roots[0]
    else:
        raise ValueError("Failed to find roots using Coppersmith's method")

#####pulled from sage implementation for small roots!
def adapted_small_roots(self, X=None, beta=1.0, epsilon=None, **kwds):
    N = self.parent().characteristic()

    if not self.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    beta = RR(beta)
    if beta <= 0.0 or beta > 1.0:
        raise ValueError("0.0 < beta <= 1.0 not satisfied.")

    f = self.change_ring(ZZ)

    P,(x,) = f.parent().objgens()

    delta = f.degree()

    if epsilon is None:
        epsilon = beta / 8

    m = max(beta**2 / (delta * epsilon), 7 * beta / delta).ceil()

    t = int((delta * m * (1 / beta - 1)).floor())

    if X is None:
        X = (0.5 * N**(beta**2 / delta - epsilon)).ceil()

    g = [x**j * N**(m-i) * f**i for i in range(m) for j in range(delta)]
    g.extend([x**i * f**m for i in range(t)])  # h

    B = Matrix(ZZ, len(g), delta * m + max(delta, t))
    for i in range(B.nrows()):
        for j in range(g[i].degree() + 1):
            B[i, j] = g[i][j] * X**j

    B = compute_and_plot_gso(B, "hastad")

    f = sum([ZZ(B[0, i] // X**i) * x**i for i in range(B.ncols())])
    R = f.roots()

    ZmodN = self.base_ring()
    roots = set([ZmodN(r) for r, m in R if abs(r) <= X])
    Nbeta = N**beta
    return [root for root in roots if N.gcd(ZZ(self(root))) >= Nbeta]

def hastads_attack_lattice(ciphertexts, moduli, e):

    try:
        recovered_message = lattice_attack(ciphertexts, moduli, e)
        return recovered_message
    except ValueError as e:
        print(e)
    