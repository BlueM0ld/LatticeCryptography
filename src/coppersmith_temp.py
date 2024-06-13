from sage.all import Integer, log, matrix, random_prime, RR, ZZ, PolynomialRing
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix

def plot_gso(log_gso_norms):
    plt.figure(figsize=(10, 6))
    for i, vec in enumerate(log_gso_norms):
        plt.plot(range(len(vec)), vec, label=f"Vector {i+1}")
    plt.ylabel("log Gram-Schmidt Norms")
    plt.title("LLL Reduction")
    plt.legend()
    plt.show()

def compute_and_plot_gso(M):
    reducedL = M.LLL()
    fpylll_matrix = convert_to_fpylll(reducedL)
    M = GSO.Mat(fpylll_matrix)
    M.update_gso()
    square_gso_norms = M.r()
    log_gso_norms = [RR(log(square_gso_norm, 2)/2) for square_gso_norm in square_gso_norms]
    
    plot_gso([log_gso_norms])
    
    return reducedL

def convert_to_fpylll(mat):
    return IntegerMatrix.from_matrix(mat)

def coppersmith_univariate(N, c, known_prefix, unknown_length, e=3):
    a = Integer(known_prefix, base=35)
    X = Integer('1' + '0' * unknown_length, base=35)
    
    M = matrix(ZZ, [
        [X**e, e*X**(e-1)*a, e*X*(a**(e-1)), a**e - c],
        [0, N*X**(e-1), 0, 0],
        [0, 0, N*X, 0],
        [0, 0, 0, N]
    ])

    B = compute_and_plot_gso(M)
    
    PolynomialRingZZ = PolynomialRing(ZZ, 'x')
    x = PolynomialRingZZ.gen()
    
    Q = B[0][0] * x**e // X**e + B[0][1] * x**(e-1) // X**(e-1) + B[0][2] * x // X + B[0][3]
    
    roots = Q.roots(ring=ZZ)
    
    if roots:
        recovered_root = roots[0][0]
        recovered_message = Integer(recovered_root).str(base=35)
        
        return Q, recovered_root, recovered_message
    else:
        return Q, None, None
