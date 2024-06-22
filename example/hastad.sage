from sage.all import *
from gmpy2 import iroot

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
# https://theory.stanford.edu/~gdurf/durfee-thesis-phd.pdf

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

    roots = g.small_roots()

    if roots:
        return roots[0]
    else:
        raise ValueError("Failed to find roots using Coppersmith's method")

if __name__ == "__main__":
    
    bits = 250
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

