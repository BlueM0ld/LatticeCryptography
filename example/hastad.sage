from sage.all import *
from gmpy2 import iroot
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix
from graph_plotting import compute_and_plot_gso, convert_to_fpylll 
import binascii

def generate_rsa_instance(bits=150, e=3):
    while True:
        p = random_prime(2** bits)
        q = random_prime(2** bits)
        N = p * q
        phi_N = (p-1) *(q-1)
        if gcd(phi_N, e) == 1:
            break
    return N, e

def encrypt_message(message, N , e):
    if N is None:
        generate_rsa_instance(150, 3)
    
    
    message = Integer(message, base=35)
    cipher_text = message ** e % N 
    return cipher_text, N

# encrypt message with linear padding
def encrypt_message_lb(message, N , e):
    m = Integer(int(message, base=35))
    random_pad_factor = randint(1, N - 1)
    random_pad_offset = randint(1, N - 1)
    padded_message = random_pad_factor * m + random_pad_offset
    cipher_text = padded_message**e % N
    return (random_pad_factor, random_pad_offset, cipher_text), N


# Hastad attack we send a message to p people.
# we encrypt message with small e and different moduli N  
# The attack states that as soon as p >= e
# we can recover the message using (crt) Chinese Remainder Theorem

def hastads_attack(ciphertexts, moduli, e):

    ciphertexts = [c[2] for c in ciphertexts]

    assert len(ciphertexts) == len(moduli), "The number of ciphertexts and moduli must be the same."
    
    # use the Chinese Remainder Theorem
    combined_ciphertext = crt(ciphertexts, moduli)
    
    # compute the cube root of the combined ciphertext
    k = len(ciphertexts)
    message, exact = iroot(Integer(combined_ciphertext), k)

    if exact:
        return int(message)
    else:
        raise ValueError("Failed to recover the message, the root is not exact.")



# Hastad attack with coppersmith externsion
def hastads_attack_lattice(ciphertexts, moduli, e):
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

def adapted_small_roots(self, X=None, beta=1.0, epsilon=None, **kwds):
    r"""
    --- adapted from sage library this is not my implementation!

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
        epsilon = beta/8
    #verbose("epsilon = %f"%epsilon, level=2)

    m = max(beta**2/(delta * epsilon), 7*beta/delta).ceil()
    #verbose("m = %d"%m, level=2)

    t = int( ( delta*m*(1/beta -1) ).floor() )
    #verbose("t = %d"%t, level=2)

    if X is None:
        X = (0.5 * N**(beta**2/delta - epsilon)).ceil()
    #verbose("X = %s"%X, level=2)

    # we could do this much faster, but this is a cheap step
    # compared to LLL
    g  = [x**j * N**(m-i) * f**i for i in range(m) for j in range(delta) ]
    g.extend([x**i * f**m for i in range(t)]) # h

    B = Matrix(ZZ, len(g), delta*m + max(delta,t) )
    for i in range(B.nrows()):
        for j in range( g[i].degree()+1 ):
            B[i,j] = g[i][j]*X**j

    #B =  B.LLL(**kwds)
    B = compute_and_plot_gso(B, "hastad")

    f = sum([ZZ(B[0,i]//X**i)*x**i for i in range(B.ncols())])
    R = f.roots()

    ZmodN = self.base_ring()
    roots = set([ZmodN(r) for r,m in R if abs(r) <= X])
    Nbeta = N**beta
    return [root for root in roots if N.gcd(ZZ(self(root))) >= Nbeta]

if __name__ == "__main__":
    
    bits = 512 ## 512 we run into a number issue # at 50 MATGSO.r() return inf
    e = 3

    rsa_instances = [generate_rsa_instance(bits, e) for _ in range(e)]
     
    ciphertexts = []
    moduli = []
    for N, e in rsa_instances:
        cipher_text, N_modulus = encrypt_message_lb("messagecheck", N, e)
        ciphertexts.append(cipher_text)
        moduli.append(N_modulus)

    print(f"e = {e}")
    for i, (cipher, N) in enumerate(zip(ciphertexts, moduli)):
        print(f"ciphertext {i}: {cipher[2]}, N={N}")

    #Run Hastad
    try:
        recovered_message = hastads_attack_lattice(ciphertexts, moduli,e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

