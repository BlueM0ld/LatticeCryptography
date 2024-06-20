from sage.all import crt, Integer
from gmpy2 import iroot

def generate_rsa_instance(bits=150, e=3):
    p = random_prime(2**bits)
    q = random_prime(2**bits)
    N = p * q

    return N, e

def encrypt_message(message, N , e):
    if N is None:
        generate_rsa_instance(150, 3)
    
    message = Integer(message, base=35)
    cipher_text = message ** e % N 
    return cipher_text, N

def hastads_attack(ciphertexts, moduli, e):

    assert len(ciphertexts) == len(moduli), "The number of ciphertexts and moduli must be the same."
    
    # Combine the ciphertexts using the Chinese Remainder Theorem
    combined_ciphertext = crt(ciphertexts, moduli)
    
    # Compute the cube root of the combined ciphertext
    k = len(ciphertexts)
    message, exact = iroot(Integer(combined_ciphertext), k)
    #recovered_message = combined_ciphertext.nth_root(e)
    if exact:
        return int(message)
    else:
        raise ValueError("Failed to recover the message, the root is not exact.")

# Hastad attack we send a message to p people.
# we encrypt message with small e and different moduli N  
# The attack states that as soon as p >= e
# we can recover the message using (crt) Chinese Remainder Theorem
# https://theory.stanford.edu/~gdurf/durfee-thesis-phd.pdf

if __name__ == "__main__":
    
    bits = 250
    e = 3

    rsa_instances = [generate_rsa_instance(bits, e) for _ in range(e)]
     
    ciphertexts = []
    moduli = []
    for N, e in rsa_instances:
        cipher_text, N_modulus = encrypt_message("messagecheck", N, e)
        ciphertexts.append(cipher_text)
        moduli.append(N_modulus)

    for i, (cipher, N) in enumerate(zip(ciphertexts, moduli)):
        print(f"ciphertext {i+1}: {cipher}, N: {N}")

    #Run Hastad
    try:
        recovered_message = hastads_attack(ciphertexts, moduli,e)
        print(f"Recovered message: {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

