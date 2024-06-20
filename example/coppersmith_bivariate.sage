from sage.all import *
from fpylll import IntegerMatrix, LLL

def coppersmith_bivariate(polynomial, modulus, m_param, t_param, X_bound, Y_bound):
    # initialize polynomial ring and quotient ring
    PR = PolynomialRing(ZZ, names=('u', 'x', 'y'))
    u, x, y = PR.gens()
    quotient_ring = PR.quotient(x * y + 1 - u)
    lifted_polynomial = quotient_ring(polynomial).lift()
    U_bound = X_bound * Y_bound + 1

    # create shifts of the polynomial
    polynomial_shifts = []
    for k in range(m_param + 1):
        for i in range(m_param - k + 1):
            x_shift = x**i * modulus**(m_param - k) * lifted_polynomial(u, x, y)**k
            polynomial_shifts.append(x_shift)
    polynomial_shifts.sort()

    # find monomials
    monomials = []
    for polynomial in polynomial_shifts:
        for monomial in polynomial.monomials():
            if monomial not in monomials:
                monomials.append(monomial)
    monomials.sort()

    # Additional polynomial shifts for y
    for j in range(1, t_param + 1):
        for k in range(floor(m_param / t_param) * j, m_param + 1):
            y_shift = y**j * lifted_polynomial(u, x, y)**k * modulus**(m_param - k)
            y_shift = quotient_ring(y_shift).lift()
            polynomial_shifts.append(y_shift)
            monomials.append(u**k * y**j)

    # create the lattice
    dimension = len(monomials)
    lattice_basis = Matrix(ZZ, dimension)
    for i in range(dimension):
        lattice_basis[i, 0] = polynomial_shifts[i](0, 0, 0)
        for j in range(1, i + 1):
            if monomials[j] in polynomial_shifts[i].monomials():
                lattice_basis[i, j] = polynomial_shifts[i].monomial_coefficient(monomials[j]) * monomials[j](U_bound, X_bound, Y_bound)

    # run LLL 
    lattice_basis = IntegerMatrix.from_matrix(lattice_basis)
    LLL.reduction(lattice_basis)

    # create polynomials from the reduced basis
    PR = PolynomialRing(ZZ, names=('x', 'y'))
    x, y = PR.gens()
    polynomials = []
    for i in range(dimension):
        polynomial = 0
        for j in range(dimension):
            polynomial += monomials[j](x * y + 1, x, y) * lattice_basis[i, j] // monomials[j](U_bound, X_bound, Y_bound)
        polynomials.append(polynomial)

    # find two polynomials with gcd 1
    polynomial1 = polynomial2 = 0
    found = False
    for i, polynomial in enumerate(polynomials):
        if found:
            break
        for j in range(i + 1, len(polynomials)):
            if gcd(polynomial, polynomials[j]) == 1:
                polynomial1 = polynomial
                polynomial2 = polynomials[j]
                found = True
                break

    if polynomial1 == polynomial2 == 0:
        return 0, 0

    # calc resultant and find roots
    PR = PolynomialRing(ZZ, names='x')
    x = PR.gen()
    resultant = polynomial1.resultant(polynomial2)
    resultant = resultant(x, x)

    roots_y = resultant.roots()
    if not roots_y:
        return 0, 0
    
    solution_y = roots_y[0][0]
    intermediate_polynomial = polynomial1(x, solution_y)
    roots_x = intermediate_polynomial.roots()
    if not roots_x:
        return 0, 0
    
    solution_x = roots_x[0][0]

    return solution_x, solution_y

def main():
    # Generate RSA parameters
    p = random_prime(2**150)
    q = random_prime(2**150)
    N = p * q

    # Choose a small private exponent d
    delta = 0.27
    d = int(N**delta)
    if d % 2 == 0:
        d += 1
    while gcd(d, (p - 1) * (q - 1)) != 1:
        d += 2

    e = inverse_mod(d, (p - 1) * (q - 1))

    print(f"Pub E: {e}")
    print(f"PrI D: {d}")
    print(f"N: {N}")

    #  polynomial for Coppersmith's method
    PR = PolynomialRing(Zmod(e), names=('x', 'y'))
    x, y = PR.gens()
    A_constant = int((N + 1) // 2)
    polynomial = 1 + x * (A_constant + y)

    # conditions for Coppersmith bivariate
    X_bound = 2 * floor(N**delta)
    Y_bound = floor(N**0.5)
    m_param = 7
    t_param = int((1 - 2 * delta) * m_param)

    # run bivariate method
    found_x, found_y = coppersmith_bivariate(polynomial, e, m_param, t_param, X_bound, Y_bound)

    print(f"Found x: {found_x}")
    print(f"Found y: {found_y}")

    # get the private key from the solutions
    if found_x and found_y:
        d_recovered = (e * found_x - 1) // (N + 1 - found_y)
        print(f"Recovered private exponent (d): {d_recovered}")

        # Verify the recovered private key
        if d == d_recovered:
            print("Successfully recovered the private key!")
        else:
            print("Failed to recover the correct private key.")
    else:
        print("Failed to find the solutions.")

if __name__ == "__main__":
    main()
