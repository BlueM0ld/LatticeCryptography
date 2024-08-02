import sys
import os

# Add the path to the src directory to sys.path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(src_path)

from sage.all import Integer, random_prime
from coppersmith_univariate import coppersmith_univariate
from graph_plotting import combine_gso_norms, generate_graphs

# RSA instance setup bits
start_size = 3
end_size = 80
inc = 2

def setup_rsa_and_cipher(known_prefix, unknown_length, e, m, bit_length):
    N = random_prime(2**bit_length/2) * random_prime(2**bit_length/2)
    print("e = ",e)
    message = Integer(m, base=35)  # 0-9 - a-z
    c = message**e % N
    return N, c, known_prefix, unknown_length, e


def experiment_coppersmith_increase_bits(red):
    print("Starting experiment - Increasing number of bits ...")
    bits = [150,256] 
    
    for nbits in bits:
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher('secr00', 2, 27, 'secret', nbits)
        try:
            recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, red)
            print(f"Recovered message with N bits {nbits}: {recovered_message}")
        except Exception as ex:
            print(f"Failed with bit size {nbits}: {ex}")
    


def experiment_coppersmith_increase_bits_det(red):
    print("Starting experiment - Increasing number of bits ...")
    start = 150
    end = 300
    inc = 10
    
    for nbits in range(start, end , inc):
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher('secr00', 2, 27, 'secret', nbits)
        try:
            recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, red)
            print(f"Recovered message with N bits {nbits}")
        except Exception as ex:
            print(f"Failed with bit size {nbits}: {ex}")
    

def experiment_coppersmith_change_known_part(red):
    print("Starting experiment - Decreasing known prefix ...")
    #message = 'secret'
    #message = 'iknowyoursecretcharlie'
    #message = 'thereisanimposteramongus'
    message = os.urandom(10).hex()

      # Lengths of the unknown part
    
    for unknown_length in range(1,len(message)-1,1):
        known_prefix = message[:-unknown_length]
        nbits = 256  # Fixed bit length
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(known_prefix, unknown_length, 3, message, nbits)
        
        try:
            recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, red)
            print(f"Recovered message with known prefix: {known_prefix}")
            print(f"Was the recovered message correct? : {recovered_message == message}")
        except Exception as ex:
            print(f"Failed with known prefix length {len(known_prefix)}: {ex}")
    

def experiment_coppersmith_change_e():
    print("Starting experiment 1 ...")
    for size in range(start_size, end_size+1,inc):
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher('secr00', 2, size, 'secret', 256)#1024, 750
        recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, "BKZ")


if __name__ == '__main__':

   #run one at a time 

   #experiment_coppersmith_increase_bits("BKZ60")
   #experiment_coppersmith_increase_bits_det("LLL")
   #experiment_coppersmith_change_known_part("BKZ60")
   #experiment_coppersmith_change_e("LLL")


   #if you want to see the graphs pass 2 params the attack and the reduction type

   #combine_gso_norms("coppersmith_univariate", "BKZ60")

   print("done")