#Coron, JS. (2007). Finding Small Roots of Bivariate Integer Polynomial Equations: A Direct Approach. In: Menezes, A. (eds) Advances in Cryptology - CRYPTO 2007. CRYPTO 2007. Lecture Notes in Computer Science, vol 4622. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-74143-5_21

from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

def coron_bivariate_direct(N,k):

    X = 30
    Y = 20
    p = PolynomialRing(ZZ, 'x, y')
    x, y = p.gens()

    pxy = 127*x*y + 1207*x + 1461*y + 21


    poly_degree = max(pxy.degrees())  # total_degree shouldnt be the sum of the degress but it is bug
    print("d",poly_degree)
    i0, j0 = 1, 1  # stub numbers

    S = Matrix(ZZ, k**2, k**2, 0)
    print("Matrix S",S)

    #sxy = (x**a)*(y**b)
    #rxy = (x**i * y**i) * det_M

    #####generate N

    #  k^2 polynomials sa,b(x, y)
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


    #there are k^2 + (k + δ)^2 such polynomials
    print("no. of polynomials", k**2+(k + poly_degree)**2)
    # Linear combinations
    # https://www.youtube.com/watch?v=YwywC7s17E8

    # h(x, y) be a linear combination of the polynomials sa,b(x, y) and ri,j (x, y)


    M_rows = k**2 + (k + poly_degree)**2 # 4 + 9 = 13
    M_columns = (k+poly_degree)**2 
    M = Matrix(ZZ, M_rows, M_columns )
   
    print("there should be (d + k)^ 2 polynomials", (poly_degree + k)**2)

    #sxy = (x**a)*(y**b)
    #rxy = (x**i * y**i) * n

    Spolynomials = []
    for a in range(k):
        for b in range(k):
            sxy = (X*x**a)*(Y*y**b)*pxy(X*x,Y*y)
            Spolynomials.append(sxy)
            #print("a",a)
            #print("b",b)
            #print("sxy",sxy)
    Spolynomials.reverse()
    #print(Spolynomials)

    Rpolynomials = []
    for i in range(k+poly_degree):
        for j in range(k+poly_degree):
            rxy = (x**i * y**j) * n
            Rpolynomials.append(rxy)

    print("this", Rpolynomials)
    Rpolynomials.reverse()

    polySystem = Spolynomials + Rpolynomials
    print(*polySystem, sep = "\n")


    # Define the order of monomials
    monomial_order = [x**2*y**2, x**2*y, x*y**2, x*y, x**2, y**2, x, y, 1]

    # Initialize an empty matrix with zeros
    # Populate the matrix with coefficients of each polynomial
    for i, poly in enumerate(polySystem):
        terms = poly.monomials()
        coeff = poly.coefficients()
        print("terms", terms)
        #for term in terms:
        for j in range(len(terms)):
            monomial = terms[j]
            print(coeff)
            #coeff = term[1]
            ## Find the index of the monomial in monomial_order
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

    L2 = L.submatrix(k ** 2, k ** 2, (k + poly_degree) ** 2 - k ** 2, (k + poly_degree) ** 2 - k ** 2)

    print(L2)
    print(L2.dimensions())

    LRED = L2.LLL()#compute_and_plot_gso(L2, "coron_bivariate") # does LLL

    print("LRED SHORTEST", LRED)

    roots = []
    monomials = [(2, 0), (0, 2), (1, 0), (0, 1), (0, 0)]
    
    
    print(monomials)
    
    for i in range(5):
        u = (sum(coeff * x**i * y**j for coeff, (i, j) in zip(LRED[i], monomials)))
        #print(u)
        print("--------------------------------------------------")
        resultant_pu = pxy.resultant( u, variable = y)
        #print("Resultant Q(x):", resultant_pu)

        gcd_coeff = gcd(resultant_pu.coefficients())
        simplified_poly = resultant_pu // gcd_coeff
        print(simplified_poly)

        x0 = simplified_poly.univariate_polynomial().roots(ring=RR, multiplicities=True)
        #x0 = simplified_poly.groebner_basis()
        print(x0)


    return None


def generate_rsa_instance(bits, e):
    p = random_prime(2**bits) ## needs to be 2**bits but cant debug well
    q = random_prime(2**bits) ## needs to be 2**bits but cant debug well
    N = p * q
    return N, e, p, q

bits = 512 
e = 3
k = 2 # assert that it is greater than 0

N, e, p, q = generate_rsa_instance(bits, e)
print('p =', p)
print('q =', q)


result = coron_bivariate_direct(N,  k)
print("Result:", result)

#Bounds described in the paper but we will set it to 1
#XY < W 1/δ