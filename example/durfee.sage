from sage.all import *


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


#a balanced RSA modulus
bits = 20  # This gives a modulus of approximately 2048 bits

# Generate two distinct large primes p and q of the specified bit length
#p = random_prime(2**bits, lbound=2**(bits-1))
#q = random_prime(2**bits, lbound=2**(bits-1))
p= 754723 
q= 771181 

# Compute the modulus N
N = p * q

# Compute Euler's totient function φ(N)
phi_N = (p - 1) * (q - 1)

# Choose a public exponent e
d = 7
# Compute the private exponent d (the modular inverse of e modulo φ(N))
e= inverse_mod(d, phi_N)


# Display the generated values
print(f"p= {p} (bit length: {p.nbits()})")
print(f"q= {q} (bit length: {q.nbits()})")

#print(f"N: {N} (bit length: {N.nbits()})")
print(f"e: {e}")
print(f"d: {d}")
delta = 0.292
bound_d = N**delta

assert d < bound_d, "method wont work"

#common case
alpha = 1


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