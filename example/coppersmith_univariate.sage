from sage.all import Integer, log, matrix, random_prime, RR, ZZ, PolynomialRing
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix
from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

def coppersmith_univariate(N, c, known_prefix, unknown_length, e=3):
    # Convert known part to integer
    a = Integer(known_prefix, base=35)
    print("a (known part as integer):", a)
    
    # Create a large X value based on the length of the unknown part
    X = Integer('1' + '0' * unknown_length, base=35)
    print("X (value for unknown part):", X)
    
    # Construct the matrix M for LLL reduction
    M = matrix(ZZ, [
        [X**e, e*X**(e-1)*a, e*X*(a**(e-1)), a**e - c],
        [0, N*X**(e-1), 0, 0],
        [0, 0, N*X, 0],
        [0, 0, 0, N]
    ])
    print("Matrix M:\n")
    print(M.str(rep_mapping=lambda x : str(x.n(digits=2))))  


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
    #
    #p  = 1159186933806439500718855209573938792641634327
    #q  = 739251066278939253402876857394291282900202471
    p = random_prime(2**150)
    q = random_prime(2**150)
    N = p * q
    print("p (prime factor of N) =", p)
    print("q (prime factor of N) =", q)
    print("N (RSA modulus) = ", N)
    
    message = Integer('theenemyishere', base=35)
    print("message_int", message)

    # Encrypt with e=3
    e = 3
    c = message**e % N
    print("Encrypted message c:", c)
    print("e:", e)


    # Define the known and unknown parts of the message
    known_prefix = 'theenemyis0000'
    unknown_length = 4  # Length of the unknown part

    Q, recovered_root, recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e)

    if recovered_root is not None:
        print("Q:", Q)
        print("Recovered root:", recovered_root)
        print("Recovered Message:", recovered_message)
    else:
        print("No roots found for the polynomial Q(x).")
