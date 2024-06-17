from sage.all import *
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix

# Graph setup
def plot_gso(log_gso_norms):
    plt.figure(figsize=(10, 6))
    for i, vec in enumerate(log_gso_norms):
        plt.plot(range(len(vec)), vec, label=f"Vector {i+1}")
    plt.ylabel("log Gram-Schmidt Norms")
    plt.title("LLL Reduction")
    plt.legend()
    plt.show()
#def compute_and_plot_gso(M):
    reducedL = M.LLL()
    fpylll_matrix = convert_to_fpylll(reducedL)
    M = GSO.Mat(fpylll_matrix)
    M.update_gso()
    square_gso_norms = M.r()
    log_gso_norms = [RR(log(square_gso_norm, 2)/2) for square_gso_norm in square_gso_norms]
    
    plot_gso([log_gso_norms])
    
    return reducedL

# convert a Sage matrix to fpylll matrix
def convert_to_fpylll(mat):
    return IntegerMatrix.from_matrix(mat)

def normalise_polynomial(polynomial):
    polynomial /= polynomial.coefficients().pop(0)
    return polynomial.change_ring(ZZ)

def generate_shifted_polynomials(polynomial, N, m, d):
    R = polynomial.base_ring()
    G = []
    
    def generate_shifts(degree, variables_count):
        shifts = []
        for i in range(degree):
            for j in range(variables_count):
                shift = [0] * variables_count
                shift[j] = i
                shifts.append(shift)
        return shifts
    
    shifts = generate_shifts(d, polynomial.nvariables())
    
    for i in range(m + 1):
        base = N ** (m - i) * polynomial ** i
        for shift in shifts:
            g = base * prod(v**s for v, s in zip(polynomial.variables(), shift))
            G.append(g)
    
    return G

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
    factors = [monomial(*bounds) for monomial in monomials]
    
    for i, factor in enumerate(factors):
        lattice_basis_mtrx.rescale_col(i, factor)
    
    lattice_basis_mtrx = lattice_basis_mtrx.dense_matrix()
    lattice_basis_mtrx = compute_and_plot_gso(lattice_basis_mtrx)
    lattice_basis_mtrx = lattice_basis_mtrx.change_ring(QQ)
    
    for i, factor in enumerate(factors):
        lattice_basis_mtrx.rescale_col(i, 1 / factor)
    
    return lattice_basis_mtrx, monomials

def extract_roots(lattice_basis_mtrx, monomials, polynomial, R):
    H = Sequence([], polynomial.parent().change_ring(QQ))
    
    for h in filter(None, lattice_basis_mtrx * monomials):
        H.append(h)
        I = H.ideal()
        
        if I.dimension() == -1:
            H.pop()
        elif I.dimension() == 0:
            roots = []
            for root in I.variety(ring=ZZ):
                root = tuple(R(root[var]) for var in polynomial.variables())
                roots.append(root)
            return roots
    
    return []

def small_roots(polynomial, bounds, m=1, d=None):
    if not d:
        d = polynomial.degree()
    
    if isinstance(polynomial, Polynomial):
        x, = polygens(polynomial.base_ring(), polynomial.variable_name(), 1)
        polynomial = polynomial(x)
    
    R = polynomial.base_ring()
    N = R.cardinality()
    
    polynomial = normalise_polynomial(polynomial)
    G = generate_shifted_polynomials(polynomial, N, m, d)
    
    lattice_basis_mtrx, monomials = construct_lattice(G)
    lattice_basis_mtrx, monomials = rescale_and_reduce_matrix(lattice_basis_mtrx, monomials, bounds)
    
    return extract_roots(lattice_basis_mtrx, monomials, polynomial, R)

def coppersmith_univariate(N, c, known_prefix, unknown_length, e):
    known_prefix_int = Integer(known_prefix, base=35)
    X = Integer(35) ** unknown_length

    P.<x> = PolynomialRing(Zmod(N), 1)
    polynomial = (known_prefix_int + x) ^ e - c

    #call small roots functions
    roots = small_roots(polynomial, (X,))
    if roots:
        recovered_root = roots[0][0]
        recovered_message = known_prefix_int + recovered_root
        recovered_message_str = Integer(recovered_message).str(base=35)
        print("Recovered message:", recovered_message_str)
        return recovered_message
    else:
        print("No roots found")
        return None
