#Coron, JS. (2007). Finding Small Roots of Bivariate Integer Polynomial Equations: A Direct Approach. In: Menezes, A. (eds) Advances in Cryptology - CRYPTO 2007. CRYPTO 2007. Lecture Notes in Computer Science, vol 4622. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-74143-5_21

from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

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

    # Reverse the order 
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
                print("Matrix S MONOMILA generated:", monomial)

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
    #assert n == 127**4, "matrix S is NOT invertible"

    return n
  

def find_roots(LRED, monomial_order, k, pxy):
    roots_x0 = []

    for r in range(LRED.nrows()):
        print("--------------------------------------------------")

        bound_order = [monomial.subs({x: X, y: Y}) for monomial in monomial_order]
       # u = sum(coeff * monomial for coeff, monomial in zip(LRED[r], monomial_order[k**2:]))
        u = sum(int(coeff / bounds)* monomial for coeff, monomial, bounds in zip(LRED[r], monomial_order[k**2:], bound_order[k**2:]))

        
        gcd_coeff = gcd(u.coefficients())
        u = u // gcd_coeff

        print("P(x):", pxy)
        print("H(x):", u)
        resultant_pu = pxy.resultant(u, variable=y)

        gcd_coeff = gcd(resultant_pu.coefficients())
        simplified_poly = resultant_pu // gcd_coeff

        x0 = simplified_poly.univariate_polynomial()
        print("Q(x):", x0)

        # Find the roots of the simplified polynomial Q(x)
        roots = x0.roots(ring=ZZ, multiplicities=False)
        if roots:
            print("Roots:", roots)
            roots_x0.extend(roots)

    roots_x0 = list(set(roots_x0))
    print("x0_roots:", roots_x0)

    sol = []
    # Extract y for each root of x0
    for root in roots_x0:
        p_y = pxy(x=root)

        # Find the roots of p_y
        roots_y = p_y.univariate_polynomial().roots(ring=ZZ, multiplicities=False)
        print(f"Roots in y for x = {root}:", roots_y)
        if roots_y:
            print(root)
            print(roots_y[0])
            sol.append([root, roots_y[0]])

    return sol


def coron_d(pxy, X, Y, k):

    print("k",k)

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
            sxy = (X*x**a)*(Y*y**b)*pxy(X*x,Y*y)
            #sxy = (x**a)*(y**b)*pxy(x,y) ### found it issue with bounds!!! need to account in the lattice which i didnt!!!
            Spolynomials.append(sxy)
    Spolynomials.reverse()

    Rpolynomials = []
    for i in range(k+poly_degree):
        for j in range(k+poly_degree):
            rxy = ((X*x)**i) * ((Y*y)**j) * n
            #rxy = (x**i) * (y**j) * n

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
    print(M.str(rep_mapping=lambda x : str(x.n(digits=2))))  


    print(M.dimensions())

    L = M
    L = L.echelon_form() ### As its over the integers echelon form is HNF(hermite normal form)!

    print("Matrix L generated:")
    print(L.str(rep_mapping=lambda x : str(x.n(digits=2))))  
 
    print(L.dimensions())

    lastpos = L.nonzero_positions_in_column(L.ncols()-1)[-1] ## for the l2 matrix

    start_row = k ** 2
    start_col = k ** 2
    num_rows = (k + poly_degree) ** 2 - k ** 2
    num_cols = (k + poly_degree) ** 2 - k ** 2

    # Extract the submatrix
    L2 = L.submatrix(start_row, start_col, num_rows, num_cols)
    print("L2--------------------------------------------------")

    print(L2.str(rep_mapping=lambda x : str(x.n(digits=2))))  
 

    LRED = compute_and_plot_gso(L2, "coron_bivariate")

    print("REDUCED--------------------------------------------------")

    print(LRED.str(rep_mapping=lambda x : str(x.n(digits=2))))  

    res = find_roots(LRED, monomial_order, k, pxy)

    return res



#Example set in gailbraiths 19.3 bivariate section

P = PolynomialRing(ZZ, 'x, y')
x, y = P.gens()
    

X = 25
Y = 10
#pxy = 127*x*y - 1207*x - 1461*y + 21
pxy = 131*x*y - 1400*x + 20*y - 1286
 

poly_degree = max(pxy.degrees())  # total_degree shouldnt be the sum of the degress but it is bug
print("degree of polynomial ",poly_degree)



# Compute W
W = 0


# Iterate through the coefficients and compute the maximum value
for i in range(poly_degree + 1):
    for j in range(poly_degree + 1):
        F = abs(pxy(X**i, Y**j))
        term_value = F* X**i * Y**j
        print(term_value)
        if term_value > W:
            W = term_value

deg = (2/(3*poly_degree))
# Print the result
print("X*Y =", X*Y)
print(f"W = {W**deg}")    
assert X*Y < W**deg, "there is not a solution "
     
res = coron_d(pxy, X, Y, k=2)

print(res)