def generate_monomial_order(k):
    #x, y = var('x y ')
    monomial_order = []

    # Add (0, 0) first
    monomial_order.append(1)

    # Add (0, i) and (i, 0) for i from 1 to k
    for i in range(1, k + 1):
        monomial_order.append(y**i)
        monomial_order.append(x**i)

    # Add (m, n) for m and n from 1 to k
    for m in range(1, k + 1):
        for n in range(1, k + 1):
            monomial_order.append(x**m * y**n)

    # Reverse the order as per the requirement
    monomial_order.reverse()

    return monomial_order

def generate_lattice(k, poly_degree):
    M_rows = k**2 + (k + poly_degree)**2 # 4 + 9 = 13
    M_columns = (k+poly_degree)**2 
    M = Matrix(ZZ, M_rows, M_columns)

    return M


def generate_n(k):

    i0 = 1
    j0 = 1
    S = Matrix(ZZ, k**2, k**2, 0)

    polynomials = []
    for a in range(k):
        for b in range(k):
            sxy = (x**a)*(y**b)*pxy
            polynomials.append(sxy)
    
    print(polynomials)

    # matrix S with coefficients of specified monomials
    for index, poly in enumerate(polynomials):
        row = []
        for i in range(k):
            for j in range(k):
                monomial = x^(i0 + i) * y^(j0 + j)
                coeff = poly.coefficient(monomial)
                row.append(coeff)
        S[index] = row

    S = S.transpose()
    print("Matrix S generated:")
    print(S.str(rep_mapping=lambda x : str(x.n(digits=2))))  
    
    # determinant of the matrix M
    n = abs(det(S)) ## n := | det S|
    
    print("|Det(s)|: ", n)
    assert n > 0, "matrix S is NOT invertible"
    assert n == 127**4, "matrix S is NOT invertible"

    return n
  


def coron_d(pxy, X, Y, k):

    poly_degree = max(pxy.degrees())  # total_degree shouldnt be the sum of the degress but it is bug
    print("degree of polynomial ",poly_degree)
    i0, j0 = 1, 1  # stub numbers

    n = generate_n(k)

    M_rows = k**2 + (k + poly_degree)**2
    M_columns = (k+poly_degree)**2 
    M = Matrix(ZZ, M_rows, M_columns )
   
    #TODO: BIGGGG NEED TO ACCOUNT FOR BOUNDS 

    Spolynomials = []
    for a in range(k):
        for b in range(k):
            #sxy = (X*x**a)*(Y*y**b)*pxy(X*x,Y*y)
            sxy = (x**a)*(y**b)*pxy(x,y) ### found it issue with bounds!!! need to account in the lattice which i didnt!!!
            Spolynomials.append(sxy)
    Spolynomials.reverse()

    Rpolynomials = []
    for i in range(k+poly_degree):
        for j in range(k+poly_degree):
            rxy = (x**i) * (y**j) * n
            Rpolynomials.append(rxy)

    print("this", Rpolynomials)
    Rpolynomials.reverse()

    polySystem = Spolynomials + Rpolynomials
    print(*polySystem, sep = "\n")

    M = generate_lattice(k, poly_degree)

    monomial_order = generate_monomial_order(k)
    print(monomial_order)

    for i, poly in enumerate(polySystem):
        terms = poly.monomials()
        coeff = poly.coefficients()
        for j in range(len(terms)):
            monomial = terms[j]
            index = monomial_order.index(monomial)
            M[i, index] = coeff[j]

    print("Matrix M generated:")
    print(M)  

    print(M.dimensions())

    L = M
    L = L.echelon_form() ### As its over the integers echelon form is HNF(hermite normal form)!

    print("Matrix L generated:")
    print(L)  
    print(L.dimensions())

    lastpos = L.nonzero_positions_in_column(L.ncols()-1)[-1] ## for the l2 matrix

    start_row = k ** 2
    start_col = k ** 2
    num_rows = (k + poly_degree) ** 2 - k ** 2
    num_cols = (k + poly_degree) ** 2 - k ** 2

    # Extract the submatrix
    L2 = L.submatrix(start_row, start_col, num_rows, num_cols)
    print("L2--------------------------------------------------")

    print(L2)  

    LRED = L2.LLL()

    print("REDUCED--------------------------------------------------")

    print(LRED.str(rep_mapping=lambda x : str(x.n(digits=2))))  

    for r in range(LRED.nrows()):
        print("--------------------------------------------------")

        u = sum(coeff * monomial for coeff, monomial in zip(LRED[r], monomial_order[k**2:]))
        #print(u)
        
        gcd_coeff = gcd(u.coefficients())
        u = u // gcd_coeff
        
        print("P(x):", pxy)
        print("H(x):", u)
        resultant_pu = pxy.resultant(u, variable = y)

        gcd_coeff = gcd(resultant_pu.coefficients())
        simplified_poly = resultant_pu // gcd_coeff

        
        x0 = simplified_poly.univariate_polynomial()
        print("Q(x):", x0)
        factors = x0.factor()
        print("Factors:", factors)
        
        # Find the roots
        roots = x0.roots()
        print("Roots:", roots)
        
        # Compute the discriminant
        discriminant = x0.discriminant()
        print("Discriminant:", discriminant)
        
        print("-" * 50)
        #print(x0.discriminant())




bits = 512 # bitlength of primes
p = random_prime(2**bits)
q = random_prime(2**bits)
N = p*q



lbits = 350  # number of lower bits of p
ln = 2**lbits

x0 = p // ln  # upper bits of p
y0 = q // ln  # upper bits of q

P_high_bits = p % ln
Q_high_bits = (N * inverse_mod(P_high_bits, ln)) % ln
assert Q_high_bits == q % ln
X = Y = 2**(bits + 1 - lbits)  # bounds on x0 and y0

P = PolynomialRing(ZZ, 'x, y')
x, y = P.gens()
    
#pxy = (ln * x + P_high_bits) * (ln * y + Q_high_bits) - N

X = 30
Y = 20
pxy = 127*x*y - 1207*x - 1461*y + 21
print(pxy)
 

res = coron_d(pxy, X, Y, k=2)

print(res)
