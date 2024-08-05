import sys
import os

# Add the path to the src directory to sys.path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../src'))
sys.path.append(src_path)

from sage.all import Integer, random_prime
from hastad import hastads_attack_lattice
from graph_plotting import generate_graphs, combine_gso_norms


def generate_rsa_instance(bits, e):
    while True:
        p = random_prime(2**bits)
        q = random_prime(2**bits)
        N = p * q
        phi_N = (p-1)*(q-1)
        if gcd(phi_N, e) == 1:
            break
    return N, e

def encrypt_message_lb(message, N, e):
    m = Integer(int(message, base=35))
    random_pad_factor = randint(1, N - 1)
    random_pad_offset = randint(1, N - 1)
    padded_message = random_pad_factor * m + random_pad_offset
    cipher_text = padded_message**e % N
    return (random_pad_factor, random_pad_offset, cipher_text), N

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


def experiment_hastad_increase_people_k():
    e = 3
    bits = 20
    message = "secret"
    for k in  [1,2,3,4,5]:  #[10, 20, 30, 40 ,50, 100, 200]:
        ciphertexts, N = set_up_hastad(message, bits, e, k)
        try:
            recovered_message = hastads_attack_lattice(ciphertexts, N, e, "BKZ60")
            print(f"recovered message with {k} recipients: {Integer(recovered_message).str(35)}")
        except ValueError as ve:
            print(f"Failed with {k} recipients: {ve}")

def experiment_hastad_increase_e():
    k = 4 #100
    bits = 20
    message = "secret"
    for e in [3, 5, 7, 11]:
        ciphertexts, N = set_up_hastad(message, bits, e, k)
        try:
            recovered_message = hastads_attack_lattice(ciphertexts, N, e, "LLL")
            print(f"recovered message with e={e}: {Integer(recovered_message).str(35)}")
        except ValueError as ve:
            print(f"Failed with e={e}: {ve}")

def experiment_hastad_increase_bits():
    e = 3
    k = 4  
    message = "secret"
    
    #bit_sizes = [64, 128, 256, 512, 1024] 
    bits = 20 

    
    for i in range(100):
        print(f"Running experiment with {bits}-bit moduli...")
        ciphertexts, N = set_up_hastad(message, bits, e, k)
        try:
            recovered_message = hastads_attack_lattice(ciphertexts, N, e, "LLL")
            print(f"Recovered message with {bits}-bit moduli: {Integer(recovered_message).str(35)}")
        except ValueError as ve:
            print(f"Failed with {bits}-bit moduli: {ve}")
    #combine_gso_norms("hastad", "LLL")

def experiment_hastad_vary_message_length():
    e = 3 #5
    k = 4 #10  
    bits = 20 #512  
    #message_lengths = [16, 32, 64, 128, 256]  
    #message_lengths = [8, 16, 32, 64, 128]  
    message_lengths = [4, 6, 16, 32, 64]  
    for length in message_lengths:
        message = "A" * length  
        print(f"Running experiment with message length {length}...")
        
        # Setup for Hastad's attack
        ciphertexts, N = set_up_hastad(message, bits, e, k)
        
        try:
            recovered_message = hastads_attack_lattice(ciphertexts, N, e, "BKZ60")# LLL, BKZ, BK40, BKZ60
            print(f"Recovered message with length {length}: {Integer(recovered_message).str(35)}")
        except ValueError as ve:
            print(f"Failed with message length {length}: {ve}")




if __name__ == '__main__':

   #experiment_hastad_increase_bits()
   #experiment_hastad_vary_message_length()
   #experiment_hastad_increase_e()
   #experiment_hastad_increase_people_k()
   combine_gso_norms("hastad", "LLL")

   print("done")

