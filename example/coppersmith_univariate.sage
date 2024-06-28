from sage.all import Integer, log, matrix, random_prime, RR, ZZ, PolynomialRing
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix
from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

def coppersmith_univariate(N, c, known_prefix, unknown_length, e=3):
    # Convert known part to integer
    a = Integer(known_prefix, base=35)
    
    # Create a large X value based on the length of the unknown part
    X = Integer('1' + '0' * unknown_length, base=35)
    
    # Construct the matrix M for LLL reduction
    M = matrix(ZZ, [
        [X**e, e*X**(e-1)*a, e*X*(a**(e-1)), a**e - c],
        [0, N*X**(e-1), 0, 0],
        [0, 0, N*X, 0],
        [0, 0, 0, N]
    ])

    # Compute and plot GSO norms, and get the reduced basis
    B = compute_and_plot_gso(M, "coppersmith_univariate")
    
    # Define the polynomial ring and variable x
    PolynomialRingZZ = PolynomialRing(ZZ, 'x')
    x = PolynomialRingZZ.gen()
    
    # Extract the first polynomial Q from the reduced basis B
    Q = B[0][0] * x**e // X**e + B[0][1] * x**(e-1) // X**(e-1) + B[0][2] * x // X + B[0][3]
    
    # Find the roots of Q(x) over the ring of integers (ZZ)
    roots = Q.roots(ring=ZZ)
    
    if roots:
        recovered_root = roots[0][0]

        # Convert the recovered root to a base-35 string to recover the message
        recovered_message = Integer(recovered_root).str(base=35)
        
        return Q, recovered_root, recovered_message
    else:
        return Q, None, None

if __name__ == '__main__':
    # RSA modulus N
    N = random_prime(2**150) * random_prime(2**150)
    message = Integer('thepasswordfortodayisswordfish', base=35)

    # Encrypt with e=3
    e = 3
    c = message**e % N

    # Define the known and unknown parts of the message
    known_prefix = 'thepasswordfortodayis000000000'
    unknown_length = 9  # Length of the unknown part

    Q, recovered_root, recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e)

    if recovered_root is not None:
        print("Q:", Q)
        print("Recovered root:", recovered_root)
        print("Recovered Message:", recovered_message)
    else:
        print("No roots found for the polynomial Q(x).")
