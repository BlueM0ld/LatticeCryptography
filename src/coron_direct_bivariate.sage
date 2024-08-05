import itertools
from sage.all import *
import matplotlib.pyplot as plt
from fpylll import GSO, IntegerMatrix
from graph_plotting import compute_and_plot_gso,generate_graphs


def check_bound_XY(p,x,y,X,Y,d):
    W = 0

    # Iterate through all the monomials of p to calculate W
    for term in p:
        coefficient, monomial = term  
        degree_x = monomial.degree(x)  
        degree_y = monomial.degree(y)  
        term_value = coefficient * (X ** degree_x) * (Y ** degree_y)
        W = max(W, term_value)


    if(X * Y > W**(2/(3*d))): 
        print("Bounds aren't met X*Y > W**(2/(3*d), if solution is not found re-adjust bounds!")
        return False
        
    #print("X * Y",X * Y)
    #print(" W**(2/(3*d))", W**(2/(3*d)))
    return True


def generate_terms(k, d, i_0, j_0, x, y):
    # Generate the product ranges and reverse them
    LTO = [(i, j) for i in range(k) for j in range(k)][::-1]
    RTO = [(i, j) for i in range(k + d) for j in range(k + d)][::-1]
    
    # Generate the terms
    terms = [x**(i + i_0) * y**(j + j_0) for i, j in LTO]
    
    # Sort the terms in reverse order
    terms.sort(reverse=True)
    
    return terms,LTO,RTO


# Builds the matrix S and calculates n = |det S|.
def calculate_determinant_s(p,x,y,LTO,terms,k, debug=False):

    print("Caculating the determinant of Matrix S ...")

    S = Matrix(ZZ, k**2, k**2)

    for r, (a, b) in enumerate(LTO):
        s_ab = x**a * y**b * p
        coeffs, mons = zip(*list(s_ab))
        s_dict = dict(zip(mons, coeffs))
        row = vector([s_dict.get(t, 0) for t in terms])
        S[r] = row

    n = abs(det(S))
    if debug:
        print("S generated")
        print(S.str(rep_mapping=lambda x : str(x.n(digits=2))))  
    
    print("n",n)

    return n

# Construction of the matrix M
def generate_matrix_M(k,d, p , x , y , LTO,RTO, n,s_terms, X, Y):
    num_rows = k**2 + (k + d)**2
    num_cols = (k + d)**2

    M = Matrix(ZZ, num_rows, num_cols)
    

    for ri, (i, j) in enumerate(LTO):
        s_ab = x**i * y**j * p
        coeffs, mons = zip(*list(s_ab))
        s_dict = dict(zip(mons, coeffs))
        row = vector([s_dict.get(t, 0) * t(x=X, y=Y) for t in s_terms])
        M[ri] = row


    for ri, (i, j) in zip(range(k**2, num_rows), RTO):
        r_ab = x**i * y**j * n
        coeffs, mons = zip(*list(r_ab))
        r_dict = {mon: coeff for mon, coeff in zip(mons, coeffs)}
        row = vector([n * t(x=X, y=Y) if t in r_dict else 0 for t in s_terms])
        M[ri] = row
    
    return M

def coron_d(p, X, Y, k, i_0 = 1, j_0 = 0,debug=True):

    P = PolynomialRing(ZZ, 'x, y')
    x, y = P.gens()

    d = max(p.degrees())

    bound_flag = check_bound_XY(p,x,y,X,Y,d)
 

    # Generate the term order for the left hand and right hand of the matrix
    terms, LTO, RTO = generate_terms(k, d, i_0, j_0, x, y)

    # Calculate the sum of x raised to the power of i for i in the range of the degree of the first term divided by 2 plus 2
    f = sum(x**i for i in range(terms[0].degree() // 2 + 2))
    f = f * f(x=y)

    # Determine the highest order term and its degree
    highest_order = terms[0]
    hdegree = max(highest_order.degree(x), highest_order.degree(y))

    #  monomials other than the k^2 in terms
    hid_terms = [t for t in list(zip(*list(f)))[1] if max(t.degree(x), t.degree(y)) <= hdegree and t not in terms]

    # Combine the original terms with the filtered terms
    s_terms = terms + hid_terms

    n = calculate_determinant_s(p,x,y,LTO,terms,k)


    M = generate_matrix_M(k,d, p , x , y , LTO, RTO,n,s_terms, X, Y)


    if debug:
        print("M generated")
        print(M.str(rep_mapping=lambda x : str(x.n(digits=2))))  


    # generate HNF of matrix M, echelon does it :) for interger matrices
    M = M.echelon_form()

    print("Matrix L generated...")
    if debug:
        print(M.str(rep_mapping=lambda x : str(x.n(digits=2))))  
        print(M.dimensions())


    print("Matrix L2 generated sliced from Matrix M ...")
    # Get the submatrix of M that contain useful polynomials which we reduce
    L = M[list(range(k**2, k**2 + len(hid_terms))), list(range(k**2, k**2 + len(hid_terms)))]
    if debug:
        print(L.str(rep_mapping=lambda x : str(x.n(digits=2))))  
        print(L.dimensions())

    L = compute_and_plot_gso(L, "coron_direct_bivariate", red)

    print("Calculating roots...")
    h = sum(coeff * term // (X**term.degree(x) * Y**term.degree(y)) for (coeff, term) in zip(L[0], hid_terms))

    #print("h", h)

    # Takes the resultant of h(x, y) and p(x, y).
    q = h.resultant(p, variable=y)

    #print("q", q)

    roots_x = q.univariate_polynomial().roots(ring=ZZ, multiplicities=False)

    roots_y = []
    
    # Run an exhaustive search
    for x_0 in roots_x:
        p_sub_x = p(x=x_0)
        p_x = p_sub_x.univariate_polynomial()
        y_0 = p_x.roots(ring=ZZ, multiplicities=False)
            
        if y_0:
            roots_y.append(y_0[0])


    if debug and len(roots_x) > 0 and len(roots_y) > 0:
        print(f"Found roots for p:  x = {roots_x}, y = {roots_y}.")

    return roots_x, roots_y
    


def coron_attack(mask,N,p,q,reduction, k):

    global red
    red = reduction

    P = PolynomialRing(ZZ, 'x, y')
    x, y = P.gens()

    p_0 = p - (p % mask)
    q_0 = (N // p_0) - ((N // p_0) % mask)

    X = mask * 2
    Y = mask * 2

    polynomial = (x + p_0) * (y + q_0) - N

    print(f"Bounds X: {X}")
    print(f"Bounds Y: {Y}")
    print(f"Polynomial generated")
    #print(f"N = {N}")

    #call small roots functions
    roots = coron_d(polynomial, X, Y, k)


    # We need to check every combination! 
    if roots:
        x_s, y_s = roots
        root_p, root_q = 0, 0

        for x_0 in x_s:
            for y_0 in y_s:
                root_p = p_0 + x_0
                root_q = q_0 + y_0


            if root_p * root_q == N:
                print(f"found p= {root_p} and q= {root_q}")
                return root_p,root_q

    else:
        print("No roots found")
        return None


def setup(bits, debug=True):
    # Generate prime numbers P and Q that are balanced
    p = random_prime(2**bits, 2**bits -1 )
    q = random_prime(2**bits, 2**bits -1 )

    # Calculate Modulus
    N = p * q
    masked_bits = floor(1/4*(log(N.nbits(),2)))
    mask = 2**masked_bits

    # Compute the bit sizes of p and q
    p_bits = p.nbits()
    q_bits = q.nbits()
    
    # Print the bit sizes of p and q
    if debug:
        print(f"Setting up Corons Bivariate Attack ...")
        print("Bit size of N:", N.nbits())
        print("p:",p)
        print("q:",q)

    return p, q, N, mask


#def experiment_coron_1():
#    p, q, N, mask = setup(10)
#    coron_attack(mask, N,p,q , "LLL", k=2)

#experiment_coron_1()
#generate_graphs("coron_d", "LLL")