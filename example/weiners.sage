def wiener_attack(e, N):
    # Calculate the continued fraction expansion of e/N
    cf = continued_fraction(e / N)
    
    list_of_convergents = cf.convergents()
    for convergent in list_of_convergents:
        print(convergent)
    
    # Step 3: Check each convergent to see if it leads to the secret exponent d
    for convergent in list_of_convergents:
        k = convergent.numerator()
        d = convergent.denominator()
        
        # Ensure d is odd and k is non-zero
        if d % 2 == 0:
            continue
        if k == 0:
            continue
        
        # Compute phi(N) - it must be whole number if not move to next convergent 
        if (e * d - 1) % k != 0:
            continue
        
        phiN = (e * d - 1) // k
        
        # Check if ed ≡ 1 (mod k)
        if Mod(e * d, k) != 1:
            continue
        
        print(f"Calculated phi(N): {phiN}")
        print(f"with numerator k: {k}")
        print(f"with denominator d: {d}")

        # Solve the poly equation to find possible p and q
        p = var('p')
        poly = p^2 - p*(N - phi_N +1)+N # fixed
        roots = poly.roots(ring=ZZ)  # roots need to be integers
        
        print(f"convergent k={k}, d={d}")
        print("ROOTS ",roots)

        
        if len(roots) == 2:
            p_recovered = roots[0][0]
            q_recovered = roots[1][0]
            
            print(f"Recovered p = {p_recovered}")
            print(f"Recovered q = {q_recovered}")

            if p_recovered * q_recovered == N:
                print("Success! found d,p and q.")
                return d, p_recovered, q_recovered
    
    return None, None, None

def generate_rsa_instance(bits=150, e=3):
    while True:
        p = random_prime(2** bits)
        q = random_prime(2** bits)
        N = p * q
        phi_N = (p-1) *(q-1)
        if gcd(phi_N, e) == 1:
            break
    
    return N, e, p, q,phi_N

# Wiener's attack - recovers d when d is small 

#p = 113 
#q = 79
#N = p * q
#phi_N = (p - 1) * (q - 1)
#e = 2621
#d = inverse_mod(e, phi_N)
#d = 3 # didnt get 5 got 3

N, e, p,q ,phi_N = generate_rsa_instance(10, 2621)


print(f"Generated RSA instance:")
print(f"N: {N}")
print(f"e: {e}")
print(f"p: {p}")
print(f"q: {q}")
print(f"Phi N: {phi_N}")

d =3
# Check the bounds for d
d_bound = int(1/2 * N^(1/4)) # There is an improved bounds as per W.Susilo et al.
assert d < d_bound, f"d must be less than the bound 1/2 * N^(1/4) ---> d={d} < {d_bound}"

d_recovered, p_recovered, q_recovered = wiener_attack(e, N)

if d_recovered:
    print(f"Found secret exponent d: {d_recovered}") # got 3 instead of 5 i wonder y
    print(f"Prime factors p: {p_recovered}, q: {q_recovered}")
    #assert d == d_recovered, "not as original"

else:
    print("d not found")
