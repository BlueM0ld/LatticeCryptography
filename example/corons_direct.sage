def coron_bivariate_direct(N,k):

    X = 1
    Y = 1
    p = PolynomialRing(ZZ, 'x, y')
    x, y = p.gens()

    #pxy = (P_high_bits + x) * (Q_high_bits + y) - N
    pxy = (1+ x) * (1+ y) - N


    poly_degree = pxy.degree()

    i0, j0 = 1, 1  # stub numbers

    S = Matrix(ZZ, k**2, k**2, 0)
    print("Matrix S",S)

    #sxy = (x**a)*(y**b)
    #rxy = (x**i * y**i) * det_M

    #####generate N

    # Create k^2 polynomials sa,b(x, y)
    polynomials = []
    for a in range(k):
        for b in range(k):
            sxy = (x**a)*(y**b)*pxy
            polynomials.append(sxy)

    # Fill the matrix M with coefficients of specified monomials
    for index, poly in enumerate(polynomials):
        row = []
        for i in range(k):
            for j in range(k):
                monomial = x^(i0 + i) * y^(j0 + j)
                coeff = poly.coefficient(monomial)
                row.append(coeff)
        S[index] = row

    print("Matrix S generated:")
    print(S.str(rep_mapping=lambda x : str(x.n(digits=2))))
    # Calculate the determinant of the matrix M

    n = det(S)
    
    print("Det(s):", n)


def generate_rsa_instance(bits, e):
    p = random_prime(2**bits) ## needs to be 2**bits but cant debug well
    q = random_prime(2**bits) ## needs to be 2**bits but cant debug well
    N = p * q
    return N, e, p, q

bits = 150 
e = 3
k = 3 # assert that it is greater than 0

N, e, p, q = generate_rsa_instance(bits, e)
print('p =', p)
print('q =', q)


result = coron_bivariate_direct(N,  k)
print("Result:", result)

#Bounds described in the paper but we will set it to 1
#XY < W 1/δ