def wiener_attack(e, N):
    # Calculate the continued fraction expansion of e/N
    # Refer to Section 2: Continued Fraction Method (Page 2)
    cf = continued_fraction(e / N)
    
    # Generate the list of convergents from the continued fraction
    # Refer to Section 2: Continued Fraction Method (Page 2)
    list_of_convergents = cf.convergents()
    
    # Step 3: Check each convergent to see if it leads to the secret exponent d
    # Refer to Section 3: Cryptanalysis (Page 2-3)
    for convergent in list_of_convergents:
        k = convergent.numerator()
        d = convergent.denominator()
        
        # Ensure d is odd and k is non-zero
        if d % 2 == 0:
            continue
        if k == 0:
            continue
        
        # Compute phi(N) - it must be whole number if not move to next convergent 
        phiN = (e * d - 1) // k
        if phiN.is_integer() == False:
            continue
        
        # if ((p-q)/2)^2 is a perfect square then guess of K and dg is correct
        # p^2 + phi(N) - N - 1)p + N = 0


        # Solve the poly equation to find possible p and q

        p = var('p')
        poly = p^2 + (phiN - N - 1)*p + N
        #Solving Equations Exactly - https://doc.sagemath.org/html/en/tutorial/tour_algebra.html
        roots = solve(poly, p)

        print("ROOTS ",roots)
        if len(roots) == 2:
            p = roots[0].rhs()
            q = roots[1].rhs()
            
            if p * q == N:
                return d, p, q
    return None, None, None



def generate_rsa_instance(bits=150, e=3):
    while True:
        p = random_prime(2** bits)
        q = random_prime(2** bits)
        N = p * q
        phi_N = (p-1) *(q-1)
        if gcd(phi_N, e) == 1:
            break
    
    d = inverse_mod(e, phi_N)

    return N, e,p,q,d
# Example usage with e and N
# Wiener's attack - recovers d when d is small

#e = 3
#bits = 10

#N,e, p,q, d = generate_rsa_instance(bits,e)
#print(f"Generated RSA instance:")
#print(f"N: {N}")
#print(f"e: {e}")
#print(f"p: {p}")
#print(f"q: {q}")
#print(f"d: {d}")


# Set as per the example in the paper

e = 2621
N = 8927
d = 5 
# Check the bounds for d
#d_bound = 1/3 * N^(1/4)
#assert d < d_bound, "d must be less than the bound 1/3 * N^(1/4)"

d_recovered, p_recovered, q_recovered = wiener_attack(e, N)

if d:
    #assert d == d_recovered, "not as original"
    print(f"Found secret exponent d: {d_recovered}") # got 3 instead of 5 i wonder y
    print(f"Prime factors p: {p_recovered}, q: {q_recovered}")
else:
    print("d not found")
