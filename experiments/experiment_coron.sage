import sys
import os

# Add the path to the src directory to sys.path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../src'))
sys.path.append(src_path)

from sage.all import Integer, random_prime
from coron_direct_bivariate import coron_attack
from graph_plotting import combine_gso_norms, generate_graphs

# RSA instance setup bits
start_size = 150
end_size = 1024
inc = 50

def setup(bits):
    p = random_prime(2**bits)
    q = random_prime(2**bits)
    N = p * q
    masked_bits = (bits//4)
    mask = 2**masked_bits

    # Compute the bit sizes of p and q
    p_bits = p.nbits()
    q_bits = q.nbits()
    
    # Print the bit sizes of p and q
    print("Bit size of N:", N.nbits())
    
    return p, q, N, mask


def experiment_coron_increase_key_size():
    print("Starting experiment - Increasing RSA Key Sizes ...")
    start_size = 150
    end_size = 1024
    inc = 50
    for bits in range(start_size, end_size + 1, inc):
        print(f"Testing with {bits}-bit RSA keys...")
        p, q, N, mask = setup(bits)
        try:
            coron_attack(mask, N, p, q, "LLL", k=2)
        except Exception as ex:
            print(f"Failed with {bits}-bit RSA keys: {ex}")
    

def experiment_coron_increase_k_LLL():
    print("Starting experiment - Increasing Number of Lattice Reductions ...")
    
    for k in range(2, 8):
        print(f"Testing with {k} lattice reductions...")
        p, q, N, mask = setup(10)  # Fixed RSA key size for this experiment
        try:
            coron_attack(mask, N, p, q, "LLL", k)
        except Exception as ex:
            print(f"Failed with k={k}: {ex}")
    
def experiment_coron_increase_k_BKZ_40():
    print("Starting experiment - Increasing Number of Lattice Reductions ...")
    
    for k in range(2, 8):
        print(f"Testing with {k} lattice reductions...")
        p, q, N, mask = setup(10)  # Fixed RSA key size for this experiment
        try:
            coron_attack(mask, N, p, q, "BKZ40", k)
        except Exception as ex:
            print(f"Failed with k={k}: {ex}")

def experiment_coron_increase_k_BKZ_60():
    print("Starting experiment - Increasing Number of Lattice Reductions ...")
    
    for k in range(2, 8):
        print(f"Testing with {k} lattice reductions...")
        p, q, N, mask = setup(10)  # Fixed RSA key size for this experiment
        try:
            coron_attack(mask, N, p, q, "BKZ60", k)
        except Exception as ex:
            print(f"Failed with k={k}: {ex}")
    

if __name__ == '__main__':


    # If results folder exists delete it 

    #experiment_coron_increase_key_size() fails

    # Uncomment all below and execute you shold find that they all run


    # If you want to run lattice instances
    #experiment_coron_increasing_k_LLL()
    #experiment_coron_increase_k_BKZ_40()
    #experiment_coron_increase_k_BKZ_60()

    # If you want to see combined graphs execute this 

    #combine_gso_norms("coron_direct_bivariate", "LLL", None)
    #combine_gso_norms("coron_direct_bivariate", "BKZ40", None)
    #combine_gso_norms("coron_direct_bivariate", "BKZ60", None)


    # If you want to see individual graphs execute this 

    #generate_graphs("coron_direct_bivariate", "LLL", None)
    #generate_graphs("coron_direct_bivariate", "BKZ40", None)
    #generate_graphs("coron_direct_bivariate", "BKZ60", None)

   print("done")