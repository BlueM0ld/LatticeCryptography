from sage.all import *
import random

def generate_monomial_order(m,n):
    #This is set for m = 3 and t =1 but need to understand the ordering better
    monomial_order = [1]

    for i in range(1, m + 1):
        monomial_order.append(x**i)
        for j in range(1, i + 1):
            monomial_order.append((x**i) * (y**j))

    for j in range(1, m+n+1):
        monomial_order.append((x**(j-1)) * (y**j))

    #print(monomial_order)
    return monomial_order


def boneh_durfee(pxy,m,t,X,Y):

    print("attempt Boneh durfee")

    # define polynomial gxy and hxy 

    
    #PP = PolynomialRing(ZZ, 'x, y')
    #x, y = PP.gens()

    #gxy = x**i * pxy**k * e **(m-k)  # x shifts
    #hxy = y**i * pxy**k * e **(m-k)  # y shifts

    # x-shifts polynomials 
    x_shifts = []
    for k in range(m + 1):
        for i in range(m - k + 1):
            gxy = ((X*x)**i) * (pxy(X*x, Y*y)**k) * (e **(m-k)) 
#            print("x^", i)
#            print("f^", k)
#            print("e^", m-k)
#            print("-----------------------")
            x_shifts.append(gxy)
    x_shifts.sort()

    #print(x_shifts)

    y_shifts = []
    for k in range(m + 1):
        for j in range(1, t+1):
            hxy = ((Y*y)**j) * (pxy(X*x, Y*y)**k) * (e **(m-k)) 
#           print("y^", j)
#           print("f^", k)
#           print("e^", m-k)
#           print("-----------------------")
            y_shifts.append(hxy)
    #y_shifts.sort()

    #print(y_shifts)
    
    w = (m+1)*(m+2)//2 + (t)*(m+1) 


    M = Matrix(ZZ, w, w)

    combine_shifts = x_shifts + y_shifts
    #print(*combine_shifts, sep = "\n")
    
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

    #print(e)

    #leave it here and gotta figure this out ??
    DetX = e**((m*(m+1)*(m+2))/3) * X**((m*(m+1)*(m+2))/3) * Y**((m*(m+1)*(m+2))/6)

    DetY = e**((t*m*(m+1))/2) * X**((t*m*(m+1))/2) * Y**((t*(m+1)*(m+t+1))/2)

    determinant_M = det(M)
    print(determinant_M)
    if determinant_M != DetX * DetY:
        print("It aint matching! DetM != DetX*DetY ")
        
    if determinant_M < e**(m*w):
        print("this condition aint met Det(L) <e^mw")


    if DetX < e**m:
        print("It aint matching! DetX <e^m ")

    if DetX != (e**20)*(X**20)*(Y**10):
        print("Something is fishy")


    L = M.LLL()
    print(L.str(rep_mapping=lambda x : str(x.n(digits=2))))  

    # Further steps to extract solution from L
    # Extracting the polynomials from the reduced basis
    basis_polynomials = []

    print(monomial_order)

    poly = (sum(coeff * monom for coeff, monom in zip(L[0], monomial_order)))
    poly = (sum(coeff * monom for coeff, monom in zip(L[1], monomial_order)))
    basis_polynomials.append(poly)



    root_list =[] 
    for g1 in basis_polynomials:
        for g2 in basis_polynomials:
            resultant = g1.resultant(g2, variable = y)
            if(resultant.is_constant()):
                print("resultant is constant", resultant)
                continue
            else:
                print("resultant added")
                x0 = resultant.univariate_polynomial()
                gcd_coeff = gcd(x0.coefficients())
                x0 = x0 // gcd_coeff
                roots = x0.roots(ring=ZZ, multiplicities=False)
                if roots:
                    print("Roots:", roots[0])
                    root_list.append(roots[0])
    
    print("Roots list:", root_list)




def generate_rsa_instance(bits=150, e=65537):
    while True:
        p = random_prime(2** bits)
        q = random_prime(2** bits)
        N = p * q
        phi_N = (p-1) *(q-1)
        if gcd(phi_N, e) == 1:
            break
    
    return N, e, p, q,phi_N

def select_random_d(N, delta):
    a = floor((3/4) * (N ** delta))
    b =floor( N ** delta)
    d = random.randint(a, b)
    return d


#a balanced RSA modulus
bits = 512   #

# Generate two distinct large primes p and q of the specified bit length
p = random_prime(2**bits, lbound=2**(bits-1))
q = random_prime(2**bits, lbound=2**(bits-1))


N, e, p,q ,phi_N = generate_rsa_instance(bits, 3)


delta = 0.05
bound_d = N**0.25
alpha = 1

#d = ZZ(select_random_d(N, delta))
#print(f"d= {d}")

#e = inverse_mod(d, phi_N)
#d = 3
d = inverse_mod(e, phi_N)

assert e * d % phi_N == 1, "e and d are not multiplicative inverses modulo phi(N)"
assert d < bound_d, "method wont work"
assert e < N**1.875, "attack is ineffective"




P = PolynomialRing(ZZ, 'x, y')
x, y = P.gens()


#as per paper
A = (N+1) //2

# k(A+S) = 1 mod (e)
pxy = x*(A + y) -1


X = floor(e**(delta))
Y = floor(e**(0.5))


m = 5
t = 1

print(f"X= {X}")
print(f"Y= {Y}")
print(f"p= {p}")
print(f"q= {q}")
print(f"N= {N}")
res = boneh_durfee(pxy,m,t,X,Y)