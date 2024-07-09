from sage.all import *
from fpylll import IntegerMatrix, LLL, GSO
from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

def generate_matrix_M1(k, p, X, Y, x, y):
    # gamma function
    def gamma(g, h):
        return (k + 1) * g + h
    
    # beta function
    def beta(i, j):
        return (k + 1)**2 + k * i + j

    # matrix dimensions
    num_rows = (k + 1)**2
    num_cols = (k + 1)**2 + k**2
    
    # set matrix M1 with zeros
    M1 = Matrix(QQ, num_rows, num_cols)

    # left-hand diagonal block
    for g in range(k + 1):
        for h in range(k + 1):
            index = gamma(g, h)
            M1[index, index] = (X**(-g) * Y**(-h))  # Ensure fractions are correctly handled
    
    # right-hand block
    for i in range(k):
        for j in range(k):
            q_ij = (x**i * y**j * p)
            for monom, coeff in q_ij.dict().items():
                g, h = monom
                row = gamma(g, h)
                col = beta(i, j)
                M1[row, col] = coeff

    # Assertion to check the size of M1
    assert M1.nrows() == (k + 1)**2, "Number of rows in M1 is incorrect"
    assert M1.ncols() == (k + 1)**2 + k**2, "Number of columns in M1 is incorrect"

    return M1

def generate_matrix_M2(M, K):
    # Matrix dimensions
    n = M.nrows()
    m = M.ncols()
    
    submatrix_cols = K**2

    # Scan the matrix from the bottom right to the top
    for j in reversed(range(m - submatrix_cols, m)):
        non_zero_positions = M.nonzero_positions_in_column(j)
        for c in reversed(range(len(non_zero_positions))):
            row_index = non_zero_positions[c]
            if -M[row_index, j] == -1:
                continue
            M.add_multiple_of_row(row_index, non_zero_positions[-1], -M[row_index, j])

    # Swap rows if necessary to make sure the leading entry is 1
    tempN = n - 1
    for i in reversed(range(m - submatrix_cols, m)):
        if M[tempN, i] != 1:
            non_zero_positions = M.nonzero_positions_in_column(i)
            if non_zero_positions:
                M.swap_rows(tempN, non_zero_positions[0])

        tempN -= 1

    # Clear fractions in the top 2K + 1 sub-rows
    def lcm(a, b):
        return abs(a * b) // gcd(a, b)

    def find_lcm_of_denominators(row):
        denominators = [frac.denominator() for frac in row if frac != 0]
        if not denominators:
            return 1
        row_lcm = denominators[0]
        for denom in denominators[1:]:
            row_lcm = lcm(row_lcm, denom)
        return row_lcm

    for row_index in range(2 * K + 1):
        row = M.row(row_index)
        row_lcm = find_lcm_of_denominators(row)
        if row_lcm != 1:
            M[row_index] = [elem * row_lcm for elem in row]

    return M

def generate_matrix_M3(M2, k):
    num_rows = 2 * k + 1
    M3 = M2[:num_rows, :]  # top (2k + 1) rows is the sublattice of M2 which is M3 as per paper?

    # Assertion to check the size of M3
    assert M3.nrows() == 2 * k + 1, "Number of rows in M3 is incorrect"
    assert M3.ncols() == M2.ncols(), "Number of columns in M3 is incorrect"

    return M3

##Used for bounds unsure if i need this FUNCTION NOT USED
def construct_diagonal_matrix_W(k, X, Y):
    dimension = (k + 1)^2
    W = matrix(QQ, dimension, dimension)
    
    for g in range(k + 1):
        for h in range(k + 1):
            index = (k + 1) * g + h
            W[index, index] = X^g * Y^h

    # Assertion to check diagonal elements
    for i in range(dimension):
        assert W[i,i] == X^(i // (k + 1)) * Y^(i % (k + 1)), "Diagonal element of W is incorrect"

    return W

def coppersmith_bivariate(N, P_high_bits, Q_high_bits, eps, k):
    print('Bivariate')
    
    P = PolynomialRing(ZZ, 'x, y')
    x, y = P.gens()
    
    pxy = (P_high_bits + x) * (Q_high_bits + y) - N
    
    print("Finding roots of P(x, y) = (P_high_bits + x) * (Q_high_bits + y) - N")
    print(pxy)
        
    X = (P_high_bits / int( N**(1/4)+eps))
    Y = (Q_high_bits / int( N**(1/4)+eps))

    # Check if method will work
    D = max((P_high_bits * Q_high_bits - N), Q_high_bits - Y, P_high_bits * Y, X * Y)
    check = int((X*Y)**(3/2)) < D

    print("This is a check bound to verify if the method will work| will this method work: ", check) 


    M1 = generate_matrix_M1(k, pxy, X, Y, x, y) 
    print("Matrix M1 generated:")
    print(M1.str(rep_mapping=lambda x : str(x.n(digits=2))))
    print("COL:", M1.ncols())
    print("ROWS:", M1.nrows())
    print("DIMENSIONS:", M1.dimensions())

    M2 = generate_matrix_M2(M1, k)
    print("M2: Transformed M1 to M2 by elementary row operations")
    print(M2.str(rep_mapping=lambda x : str(x.n(digits=2))))
    print("COL:", M2.ncols())
    print("ROWS:", M2.nrows())
    print("DIMENSIONS:", M2.dimensions())

    M3 = generate_matrix_M3(M2, k) 

    print("Matrix M3 (Top 2k + 1 rows of M2):")
    print(M3.str(rep_mapping=lambda x : str(x.n(digits=2))))
    print("COL:", M3.ncols())
    print("ROWS:", M3.nrows())
    print("DIMENSIONS:", M3.dimensions())

    L = M3[:2 * k + 1, :2 * k + 1]
    print("Matrix L (Left-hand (2k + 1) x (2k + 1) submatrix of M3):")
    print(L.str(rep_mapping=lambda x : str(x.n(digits=2))))
    
    # Assertions to check the size of L
    # Assertions to check the size of L
    assert L.nrows() == 2 * k + 1, "Number of rows in L is incorrect"
    assert L.ncols() == 2 * k + 1, "Number of columns in L is incorrect"




    #check that size of determinant is grete that N^k/4
    #det_L = abs(det(L))
    #assert det_L > N**(k/4), "|det(L)| is not greater than N^(k/4)"

    
    #L = L.LLL()
    L = L.change_ring(ZZ)
    L = compute_and_plot_gso(L, "coppersmith_bivariate")


    roots = []
    monomials = [(i, j) for i in range(k + 1) for j in range(k + 1)]

    for i in range(L.nrows()):
        u = (sum(coeff * x**i * y**j for coeff, (i, j) in zip(L[i], monomials)))
        r = pxy.resultant(u, y)
        if r.is_constant():
            print("FOUND CONSTANT")
            print(r)
            continue
        else:
            print("FOUND SOMEHTING")
            print(r)
            root = r.univariate_polynomial().roots()
            print("ROOTS?")
            print(root)
            if root:
                roots.append(root)

    return roots

def generate_rsa_instance(bits, e):
    p = random_prime(bits) ## needs to be 2**bits but cant debug well
    q = random_prime(bits) ## needs to be 2**bits but cant debug well
    N = p * q
    return N, e, p, q

bits = 150 
e = 3
eps = 1/10
k = 3
print("K",k)
N, e, p, q = generate_rsa_instance(bits, e)
log_N = log(N,2)

high_bits_length = int(((1/4 + eps) * log_N)* p.nbits())
print("high bits length",high_bits_length)

P_high_bits = p >> (p.nbits() - high_bits_length)
Q_high_bits = q >> (q.nbits() - high_bits_length)

print("HIGH BITS p",P_high_bits)
print("HIGH BITS q",Q_high_bits)


result = coppersmith_bivariate(N, P_high_bits, Q_high_bits, eps, k)
print("Result:", result)
