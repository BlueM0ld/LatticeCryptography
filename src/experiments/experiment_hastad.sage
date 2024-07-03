import sys
import os

# Add the path to the src directory to sys.path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(src_path)

from sage.all import Integer, random_prime
from hastad import hastads_attack_lattice

def generate_rsa_instance(bits=150, e=3):
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

def set_up_hastad(message, bits, e):
    rsa_instances = [generate_rsa_instance(bits, e) for _ in range(e)]
    ciphertexts = []
    moduli = []
    for N, e in rsa_instances:
        cipher_text, N_modulus = encrypt_message_lb(message, N, e)
        ciphertexts.append(cipher_text)
        moduli.append(N_modulus)

    print(f"e = {e}")
    for i, (cipher, N) in enumerate(zip(ciphertexts, moduli)):
        print(f"ciphertext {i}: {cipher[2]}, N={N}")
    return ciphertexts, moduli

def experiment_hastad_1():
    bits = 150
    e = 3
    message = "alicedidit"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_2():
    bits = 256
    e = 3
    message = "alicedidit"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_3():
    bits = 512
    e = 3
    message = "alicedidit"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_4():
    bits = 150
    e = 3
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_5():
    #WILDCARE
    bits = 256
    e = 7
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_6():
    bits = 256
    e = 3
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_7():
    bits = 512
    e = 3
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

################# FAIL HERE HMM TypeError: unable to convert '43.5849625007212+2.26618007091360*I' to a real number
def experiment_hastad_8():
    bits = 150
    e = 5
    message = "alicedidit"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_9():
    bits = 256
    e = 5
    message = "alicedidit"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_10():
    bits = 512
    e = 5
    message = "alicedidit"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_11():
    bits = 150
    e = 5
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_12():
    bits = 256
    e = 5
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)

def experiment_hastad_13():
    bits = 512
    e = 5
    message = "theenemyisonthenorthside"
    ciphertexts, N = set_up_hastad(message, bits, e)

    try:
        recovered_message = hastads_attack_lattice(ciphertexts, N, e)
        print(f"recovered message ------------> {Integer(recovered_message).str(35)}")
    except ValueError as e:
        print(e)


if __name__ == '__main__':
    experiment_hastad_1()
    experiment_hastad_2()
    experiment_hastad_3()
    experiment_hastad_4()
    experiment_hastad_5()
    experiment_hastad_6()
    experiment_hastad_7()
    #experiment_hastad_8()
    #experiment_hastad_9()
    #experiment_hastad_10()
    #experiment_hastad_11()
    #experiment_hastad_12()
    #experiment_hastad_13()
