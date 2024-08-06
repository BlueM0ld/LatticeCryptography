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


# Sage has an implementation using small roots we open this to update the call from LLL to our implemenation
def adapted_small_roots(self, X=None, beta=1.0, epsilon=None, **kwds):
    r"""
    
    adapted from sage library this is not my implementation!

    REFERENCES:

    Don Coppersmith. *Finding a small root of a univariate modular equation.*
    In Advances in Cryptology, EuroCrypt 1996, volume 1070 of Lecture
    Notes in Computer Science, p. 155--165. Springer, 1996.
    http://cr.yp.to/bib/2001/coppersmith.pdf

    Alexander May. *New RSA Vulnerabilities Using Lattice Reduction Methods.*
    PhD thesis, University of Paderborn, 2003.
    http://www.cs.uni-paderborn.de/uploads/tx_sibibtex/bp.pdf
    """

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

    # This is calling the reduction
    B =  compute_and_plot_gso(B, "hastad", reduction=red)

    f = sum([ZZ(B[0, i] // X**i) * x**i for i in range(B.ncols())])
    R = f.roots()

    ZmodN = self.base_ring()
    roots = set([ZmodN(r) for r, m in R if abs(r) <= X])
    Nbeta = N**beta
    return [root for root in roots if N.gcd(ZZ(self(root))) >= Nbeta]

def hastads_attack_lattice(ciphertexts, moduli, e, reduction):

    # passes the reduction requested
    global red
    red = reduction
    try:
        recovered_message = lattice_attack(ciphertexts, moduli, e)
        conv_recovered_message = Integer(recovered_message).str(35)
        return conv_recovered_message
    except ValueError as e:
        print(e)


# Naiive RSA key gen
def generate_rsa_instance(bits, e):
    while True:
        p = random_prime(2**bits)
        q = random_prime(2**bits)
        N = p * q
        phi_N = (p-1)*(q-1)
        if gcd(phi_N, e) == 1:
            break
    return N, e

# Encrypt plaintext with linear padding
def encrypt_message_lb(message, N, e):
    m = Integer(int(message, base=35))
    random_pad_factor = randint(1, N - 1)
    random_pad_offset = randint(1, N - 1)
    padded_message = random_pad_factor * m + random_pad_offset
    cipher_text = padded_message**e % N
    return (random_pad_factor, random_pad_offset, cipher_text), N


# Generate multiple ciphertexts with same e but different N
def set_up_hastad(message, bits, e, k, debug=False):
    rsa_instances = [generate_rsa_instance(bits, e) for _ in range(k)]
    ciphertexts = []
    moduli = []
    for N, e in rsa_instances:
        cipher_text, N_modulus = encrypt_message_lb(message, N, e)
        ciphertexts.append(cipher_text)
        moduli.append(N_modulus)

    if debug:
        print(f"e = {e}")
        for i, (cipher, N) in enumerate(zip(ciphertexts, moduli)):
            print(f"ciphertext {i}: {cipher[2]}, N={N}")
    return ciphertexts, moduli