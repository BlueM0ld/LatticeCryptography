import sys
import os

# Add the path to the src directory to sys.path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../src'))
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

    message = os.urandom(10).hex()
    print(f"Original message: {message}")

    # Remove the last 2 characters
    modified_message = message[:-2]
    print(f"Modified message: {modified_message}")
        
    for nbits in bits:
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(modified_message, 2, 3, message, nbits)
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
    message = os.urandom(6).hex()
    #message = 'secret'

    print(f"Original message: {message}")

    # Remove the last 2 characters
    m_slice = message[:-2]
    print(f"Known message: {m_slice}")
    
    for nbits in range(start, end , inc):
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher(m_slice, 2, 27, message, nbits)
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

# Read me please ! Each block reflects a experiment instance 
# Uncomment the block and run an experiment across reduction methods
# Remember to delete the results file between each run 
# You can safely uncomment the graph code stick with combine_gso, smaller amount of graphs

#    experiment_coppersmith_increase_bits("LLL")
#    experiment_coppersmith_increase_bits("BKZ40")
#    experiment_coppersmith_increase_bits("BKZ60")



#    experiment_coppersmith_increase_bits_det("LLL")
#    experiment_coppersmith_increase_bits_det("BKZ40")
#    experiment_coppersmith_increase_bits_det("BKZ60")

#    experiment_coppersmith_change_known_part("LLL")
#    experiment_coppersmith_change_known_part("BKZ40")
#    experiment_coppersmith_change_known_part("BKZ60 ")

#    experiment_coppersmith_change_e("LLL")
#    experiment_coppersmith_change_e("BKZ40")
#    experiment_coppersmith_change_e("BKZ60")


# GRAPH GENERATION ------------------------------------------------------------------------

# If you want to see combined graphs execute this 

#    combine_gso_norms("coppersmith_univariate", "LLL", None)
#    combine_gso_norms("coppersmith_univariate", "BKZ40", None)
#    combine_gso_norms("coppersmith_univariate", "BKZ60", None)


# If you want to see individual graphs execute this 

#    generate_graphs("coppersmith_univariate", "LLL", None)
#    generate_graphs("coppersmith_univariate", "BKZ40", None)
#    generate_graphs("coppersmith_univariate", "BKZ60", None)

    print("done")