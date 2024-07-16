from sage.all import *

def generate_monomial_order(m,n):
    #This is set for m = 3 and t =1 but need to understand the ordering better
    monomial_order = [1]

    for i in range(1, m + 1):
        monomial_order.append(x**i)
        for j in range(1, i + 1):
            monomial_order.append((x**i) * (y**j))

    for j in range(1, m+n+1):
        monomial_order.append((x**(j-1)) * (y**j))

    print(monomial_order)
    return monomial_order


def boneh_durfee(pxy,m,t,X,Y):

    print("attempt Boneh durfee")

    # define polynomial gxy and hxy 

    
    PP = PolynomialRing(ZZ, 'x, y')
    x, y = PP.gens()

    #gxy = x**i * pxy**k * e **(m-k)  # x shifts
    #hxy = y**i * pxy**k * e **(m-k)  # y shifts

    # x-shifts polynomials 
    x_shifts = []
    for k in range(m + 1):
        for i in range(m - k + 1):
            gxy = (x**i) * (pxy**k) * (e **(m-k)) 
            x_shifts.append(gxy)
    x_shifts.sort()

    #print(x_shifts)

    y_shifts = []
    for k in range(m + 1):
        for j in range(1, t+1):
            hxy = (y**j) * (pxy**k) * (e **(m-k)) 
            y_shifts.append(hxy)
    y_shifts.sort()

    #print(y_shifts)
    
    #Hard code values 14 by 14 until i figure out how to set dimensions
    M = Matrix(ZZ, 14, 14 )
    combine_shifts = x_shifts + y_shifts
    print(*combine_shifts, sep = "\n")
    
    monomial_order = generate_monomial_order(k,t)

    for i, poly in enumerate(combine_shifts):
        terms = poly.monomials()
        coeff = poly.coefficients()
        for j in range(len(terms)):
            monomial = terms[j]
            index = monomial_order.index(monomial)
            M[i, index] = coeff[j]

    M.swap_rows(5,6) # This is a fixed swap but will need to account for the sorting on monomials generally
    print(M.str(rep_mapping=lambda x : str(x.n(digits=2))))  

    print(e)

    #leave it here and gotta figure this out ??
    DetX = e**((m*(m+1)*(m+2))/3) * X**((m*(m+1)*(m+2))/3) * Y**((m*(m+1)*(m+2))/6)

    DetY = e**((t*m*(m+1))/2) * X**((t*m*(m+1))/2) * Y**((t*(m+1)*(m+t+1))/2)

    determinant_M = det(M)
    print(determinant_M)
    #assert determinant_M == DetX * DetY, "assertion on L"

    L = M.LLL()
    print(L.str(rep_mapping=lambda x : str(x.n(digits=2))))  

    # Further steps to extract solution from L
    # Extracting the polynomials from the reduced basis
    basis_polynomials = []

    print(monomial_order)
    for i in range(L.nrows()):
        poly = (sum(coeff * monom for coeff, monom in zip(L[0], monomial_order)))
        basis_polynomials.append(poly)

 #   for poly in basis_polynomials:
 #       print(poly)

    PR.<q> = PolynomialRing(ZZ)

    resultants = []
    for g1 in basis_polynomials:
        for g2 in basis_polynomials:
            resultant = g1.resultant(g2, variable = y)
            if(resultant.is_constant()):
                continue
            else:
                print(resultant)
                resultants.append(resultant)

    if not resultant:
        print("There's no resultants for g1 and g2")



#a balanced RSA modulus
bits = 20  # This gives a modulus of approximately 2048 bits

# Generate two distinct large primes p and q of the specified bit length
#p = random_prime(2**bits, lbound=2**(bits-1))
#q = random_prime(2**bits, lbound=2**(bits-1))
p= 754723 
q= 771181 

# Compute the modulus N
N = p * q

phi_N = (p - 1) * (q - 1)

# Choose a public exponent e
d = 7
# Compute the private exponent d (the modular inverse of e modulo φ(N))
#e= inverse_mod(d, phi_N)

#common case
alpha = 1

e = N**alpha

# Display the generated values
print(f"p= {p} (bit length: {p.nbits()})")
print(f"q= {q} (bit length: {q.nbits()})")

#print(f"N: {N} (bit length: {N.nbits()})")
print(f"e: {e}")
print(f"d: {d}")
delta = 0.292
bound_d = N**delta

assert d < bound_d, "method wont work"


P = PolynomialRing(ZZ, 'x, y')
x, y = P.gens()


#as per paper
A = (N+1) //2

# k(A+S) = 1 mod (e)
pxy = x*(A + y) -1

#paper m = 3 and t =1 but tunable

m = 3
t = 1

#Bounds
X = e**(delta)
Y = e**(0.5)
res = boneh_durfee(pxy,m,t,X,Y)