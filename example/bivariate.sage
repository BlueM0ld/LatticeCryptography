from sage.all import *

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
    M1 = Matrix(QQ, num_rows, num_cols)
    
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

def coppersmith_bivariate(N, P_high_bits, Q_high_bits, high_bits_length):
    print('Bivariate')
    
    P.<x, y> = PolynomialRing(ZZ)
    
    f1 = (P_high_bits + x) * (Q_high_bits + y) - N
    
    print("Finding roots of f1(x, y) = (P_high_bits + x) * (Q_high_bits + y) - N")
    
    k = 3  # As per paper?
    
    # choose bounds for X and Y 
    X = 2^(high_bits_length)
    Y = 2^(high_bits_length)
    
    M = generate_matrix_M1(k, f1, X, Y, x, y)
    
    # confirm matrix
    print("Matrix M1 generated:")
    print(M.str(rep_mapping=lambda x : str(x.n(digits=3))))
    print("COL:", M.ncols())
    print("ROWS:", M.nrows())
    print("DIMENSIONS:", M.dimensions())


def generate_rsa_instance(bits=512, e=3):
    p = random_prime(2 ** (bits // 2))
    q = random_prime(2 ** (bits // 2))
    N = p * q
    return N, e, p, q

# Example RSA instance
bits = 256
e = 3
eps = 1/20
k = 3

log_eps = math.log(eps)

N, e, p, q = generate_rsa_instance(bits, e)

high_bits_length = int((1/4 + log_eps) * bits)  # Number of high bits known
P_high_bits = p >> (p.nbits() - high_bits_length)
Q_high_bits = q >> (q.nbits() - high_bits_length)

result = coppersmith_bivariate(N, P_high_bits, Q_high_bits, high_bits_length)
