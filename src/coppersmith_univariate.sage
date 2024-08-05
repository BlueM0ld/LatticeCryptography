from sage.all import *
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix
from graph_plotting import compute_and_plot_gso, combine_gso_norms


def generate_shifted_polynomials(polynomial, N, m, d):
    G = []

    n_var = polynomial.nvariables()
    shifts = [[i if j == k else 0 for k in range(n_var)] for i in range(d) for j in range(n_var)]

    for i in range(m + 1):
        base = N ** (m - i) * polynomial ** i
        for shift in shifts:
            g = base * prod(v**s for v, s in zip(polynomial.variables(), shift))
            G.append(g)
    
    return G


# Construct the lattice according to the monomial set
def construct_lattice(G):
    coefficients = []
    monomials = set()
    for g in G:
        monomials.update(g.monomials())
    
    monomials = sorted(monomials, key=lambda m: (m.degree(), m))
    
    for g in G:
        row = [g.monomial_coefficient(m) for m in monomials]
        coefficients.append(row)
    
    lattice_basis_mtrx = matrix(ZZ, coefficients)
    return lattice_basis_mtrx, vector(monomials)

def rescale_and_reduce_matrix(lattice_basis_mtrx, monomials, bounds):

    #introduce bounds 
    coeffs = []
    for monomial in monomials:
      coeffs.append(monomial(*bounds))

    for i, coeff in zip(range(len(coeffs)), coeffs):
        lattice_basis_mtrx.rescale_col(i, coeff)
    

    # Call LLL and compute the GSO norms 
    reduced_lattice_basis_m = compute_and_plot_gso(lattice_basis_mtrx, "coppersmith_univariate", reduction=red)

    reduced_lattice_basis_m = reduced_lattice_basis_m.change_ring(QQ)
    
    # Remove bounds
    for i, coeff in zip(range(len(coeffs)), coeffs):
        reduced_lattice_basis_m.rescale_col(i, coeff**(-1))
    
    return reduced_lattice_basis_m, monomials



# Use a grobner method to extract roots - we look for zero ideals as they have finite roots
def extract_roots(lattice_basis_matrix, monomials, polynomial, transform_root):
    
    # Reintroduce the lattice basis
    sequence_H = Sequence([], polynomial.parent().change_ring(QQ))
    
    products = lattice_basis_matrix * monomials
    
    for monomial_product in filter(None, products):
        sequence_H.append(monomial_product)
        
        ideal_I = sequence_H.ideal()
        
        dimension = ideal_I.dimension()
        if not dimension == 0: 
            sequence_H.pop()
        elif dimension == 0:
            roots = find_roots(ideal_I, polynomial, transform_root) 
            return roots
    
    return []

def find_roots(ideal, polynomial, transform_root):
    roots = []
    variables = polynomial.variables()
    
    for o_root in ideal.variety(ring=ZZ):
        transformed_root = tuple(transform_root(o_root[var]) for var in variables)
        roots.append(transformed_root)
    
    return roots

def small_roots(polynomial, bounds):
    d = polynomial.degree()
    
    R = polynomial.base_ring()
    N = R.cardinality()
    
    polynomial = polynomial.change_ring(ZZ)
    G = generate_shifted_polynomials(polynomial, N, 1, d)
    
    #construct the lattice
    lattice_basis_mtrx, monomials = construct_lattice(G)

    # Call LLL
    lattice_basis_mtrx, monomials = rescale_and_reduce_matrix(lattice_basis_mtrx, monomials, bounds)
    
    # Solve the reduced lattice basis 
    return extract_roots(lattice_basis_mtrx, monomials, polynomial, R)



def coppersmith_univariate(N, c, known_prefix, unknown_length, e, reduction):
    # reduction type 
    global red
    red = reduction

    # convert to an integer to be used in a polynomial setup
    known_prefix_int = Integer(known_prefix, base=35)

    # set minimum bound we can set it higher but the that makes it harder 
    X = Integer(35) ** unknown_length


    # set up coopersmith univariate instance, convert based on known and unknown part

    P.<x> = PolynomialRing(Zmod(N), 1)
    polynomial = (known_prefix_int + x) ^ e - c

    #call small roots functions
    roots = small_roots(polynomial, (X,))
    if roots:
        recovered_root = roots[0][0]
        recovered_message = known_prefix_int + recovered_root
        recovered_message_str = Integer(recovered_message).str(base=35)
        print("Recovered message:", recovered_message_str)
        return recovered_message_str
    else:
        print("No roots found")
        return None


# Generate relaxed RSA instance and ciphertext
def setup_rsa_and_cipher(known_prefix, unknown_length, e, m, bit_length):
    p = random_prime(2**bit_length)
    q = random_prime(2**bit_length)

    N = p*q
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"Number of bits of N are: {N.nbits()}")
    print("e = ",e)
    message = Integer(m, base=35)  # 0-9 - a-z
    c = message**e % N
    return N, c, known_prefix, unknown_length, e

def experiment_coppersmith_1():
    print("Starting experiment 1 ...")
    #e = 19 bit = 220
    #32,64,128,256, xxx512,1024,2048
    start_size = 53
    end_size = 65
    inc = 10
    for size in range(start_size, end_size+1,inc):
        N, c, known_prefix, unknown_length, e = setup_rsa_and_cipher('itsasecre0', 1,3, 'itsasecret', 275)
        recovered_message = coppersmith_univariate(N, c, known_prefix, unknown_length, e, "LLL")
    combine_gso_norms("coppersmith_univariate", "LLL", None)

#experiment_coppersmith_1()