from sage.all import *
from fpylll import IntegerMatrix, LLL, GSO
from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

def generate_matrix_M1(k, p, X, Y, x, y):
    # gamma function
    def gamma(g, h):
        return (k + 1) * g + h
    
    # beta function
    def beta(i, j):
        return (k + 1)^2 + k * i + j

    # matrix dimensions
    num_rows = (k + 1)^2
    num_cols = (k + 1)^2 + k^2

    # set matrix M1 with zeros
    M1 = Matrix(ZZ, num_rows, num_cols)

    # left-hand diagonal block
    for g in range(k + 1):
        for h in range(k + 1):
            index = gamma(g, h)
            M1[index, index] = X^(-g) * Y^(-h)
    
    # right-hand block
    for i in range(k):
        for j in range(k):
            q_ij = (x^i * y^j * p)
            for monom, coeff in q_ij.dict().items():
                g, h = monom
                row = gamma(g, h)
                col = beta(i, j)
                M1[row, col] = coeff

    return M1

def generate_matrix_M2(M):
    M2 = M.rref()  # Reduced Row Echelon Form
    return M2

def generate_matrix_M3(M2, k):
    num_rows = 2 * k + 1
    M3 = M2[:num_rows, :]  # top (2k + 1) rows is the sublattice of M2 which is M3 as per paper?
    return M3

def coppersmith_bivariate(N, P_high_bits, Q_high_bits, high_bits_length):
    print('Bivariate')
    
    P.<x, y> = PolynomialRing(ZZ)
    
    f1 = (P_high_bits + x) * (Q_high_bits + y) - N
    
    print("Finding roots of f1(x, y) = (P_high_bits + x) * (Q_high_bits + y) - N")
    
    k = 3  # As per paper?
    
    # choose bounds for X and Y 
    X = 2^(high_bits_length)
    Y = 2^(high_bits_length)
    
    M1= generate_matrix_M1(k, f1, X, Y, x, y)
    M2 = generate_matrix_M2(M1)
    M3 = generate_matrix_M3(M2, k)
    # confirm matrix
    print("Matrix M1 generated:")
    print(M1.str(rep_mapping=lambda x : str(x.n(digits=2))))
    print("COL:", M1.ncols())
    print("ROWS:", M1.nrows())
    print("DIMENSIONS:", M1.dimensions())
    
    print("Transformed M1 -> M2 (Reduced Row Echelon Form):")
    print(M2.str(rep_mapping=lambda x : str(x.n(digits=2))))
    print("COL:", M2.ncols())
    print("ROWS:", M2.nrows())
    print("DIMENSIONS:", M2.dimensions())

    print("Matrix M3 (Top 2k + 1 rows of M2):")
    print(M3.str(rep_mapping=lambda x : str(x.n(digits=2))))
    print("COL:", M3.ncols())
    print("ROWS:", M3.nrows())
    print("DIMENSIONS:", M3.dimensions())


    num_rows = (k + 1)^2
    W = diagonal_matrix([X^g * Y^h for g in range(k + 1) for h in range(k + 1)])
    
    WM1 = W * M1
    #print("Matrix W * M1:")
    #print(WM1.str(rep_mapping=lambda x : str(x.n(digits=2))))
    
    L = M3[:, :2 * k + 1]
    print("Matrix L (Left-hand (2k + 1) x (2k + 1) submatrix of M3):")
    print(L.str(rep_mapping=lambda x : str(x.n(digits=2))))
    
    LLL_reduction = compute_and_plot_gso(M1, "bivariate")

    # LLL_reduction = L.LLL()
    print("LLL-reduced basis of matrix L:")
    print(LLL_reduction.str(rep_mapping=lambda x : str(x.n(digits=2))))
    
    short_vector = LLL_reduction[0]
    print("Short vector in the lattice:")
    print(short_vector)
    
    u = sum(coeff * x^i * y^j for coeff, (i, j) in zip(short_vector, [(i, j) for i in range(k+1) for j in range(k+1)]))
    print("Polynomial u(z0, y0) = 0:")
    print(u)
    
    # Ensure both u and f1 are in the same ring
    R = PolynomialRing(ZZ, 'x, y')
    x, y = R.gens()
    u = R(u)
    f1 = R(f1)
    
    # Compute the resultant
    resultant_y = f1.resultant(u, y)
    
    #TODO: Find roots of the resultant polynomial roots method in sage in for univariate not bivariate hmm
    return None

def generate_rsa_instance(bits=512, e=3):
    p = random_prime(2 ** (bits // 2))
    q = random_prime(2 ** (bits // 2))
    N = p * q
    return N, e, p, q

# Example RSA instance
bits = 150
e = 3
eps = 1/20 #  (eps>0): eps is a frequently used mathematical symbol to denote an arbitrary small but still positive and non zero quantity - EP
k = 3

log_eps = math.log(eps)

N, e, p, q = generate_rsa_instance(bits, e)

high_bits_length = int((1/4 + log_eps) * bits)  # Number of high bits known
P_high_bits = p >> (p.nbits() - high_bits_length)
Q_high_bits = q >> (q.nbits() - high_bits_length)

result = coppersmith_bivariate(N, P_high_bits, Q_high_bits, high_bits_length)
print("Result:", result)
