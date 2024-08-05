from graph_plotting import compute_and_plot_gso, convert_to_fpylll 

def wiener_attack_lattice(e, N):

    X = int((1/3) * N**(1/4))
    U = int((3) * N**(3/4))
    
    L = matrix(ZZ, [[e * X, 0],[N*X, U]])
    
    print(f"Lattice basis matrix L:\n{L}")
    
    L_red = L.LLL()
    print(f"Reduced lattice basis matrix L_red:\n{L_red}")
    
    # The smallest vector in the reduced basis is the first row of L_red
    u0, x0 = L_red[0]
    print(f"u0: {u0}, x0: {x0}")
    

    #det(L) <= N**nm
    if gcd(u0, x0) != 1:
        print("The gcd of u0 and x0 is not 1, the result may not be correct.")
    
    d = u0 * inverse_mod(e, x0) % N
    print(f"Recovered private exponent d: {d}")
    
    if pow(e * d % N, 1, N) != 1:
        print("Recovered private exponent d is incorrect.")
    
    if u0**2 + x0**2 > N:
        print("Howgrave-Graham's condition is not satisfied. Retry with larger bounds.")
    
    return d


def wiener_attack(e, N):
    # Calculate the continued fraction expansion of e/N
    cf = continued_fraction(e / N)
    print(cf)
    list_of_convergents = cf.convergents()
    for convergent in list_of_convergents:
        print(convergent)
    
    list_res =[]
    # Check each convergent to see if it leads to the secret exponent d
    for convergent in list_of_convergents:
        k = convergent.numerator()
        d = convergent.denominator()
        
        # Ensure d is odd and k is non-zero
        #if d % 2 == 0:
        #    print(f"skipping convergent: {convergent}")
        #    continue
        if k == 0:
            print(f"skipping convergent: {convergent}")
            continue
        
        # Compute phi(N) - it must be whole number if not move to next convergent 
        #if (e * d - 1) % k != 0:
        #    print(f"skipping convergent: {convergent}")
        #    continue
        
        g= (e*d)%k
        print("EDG--------------------",e*d)
        print("G--------------------",g)
        print("K--------------------",k)

        phiN = ((e * d) - 1) // k
        
        # Check if ed ≡ 1 (mod k)
        if Mod(e * d, k) != 1:
            print(f"skipping convergent: {convergent}")
            continue
        
        print(f"Calculated phi(N): {phiN}")
        print(f"with numerator k: {k}")
        print(f"with denominator d: {d}")

        # Solve the poly equation to find possible p and q
        p = var('p')
        poly = p^2 - p*((N - phiN) +1)+N # fixed
        print("p ",poly)

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
                d_im = inverse_mod(e, (p_recovered-1 )*(q_recovered-1))
                print("d recovered",d_im)

                list_res.append([p_recovered,q_recovered,d])
    
    return list_res

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

e=127843
N = 192649
#N, _, p,q ,phi_N = generate_rsa_instance(10, e)

print(f"Generated RSA instance:")
print(f"N= {N}")
print(f"e= {e}")
#print(f"p= {p}")
#print(f"q= {q}")
#print(f"phi_N: {phi_N}")
#d = inverse_mod(e, phi_N)
#print(f"d: {d_in}")

d =3
# Check the bounds for d
d_bound = int(1/2 * N^(1/4)) # There is an improved bounds as per W.Susilo et al.
assert d < d_bound, f"d must be less than the bound 1/2 * N^(1/4) ---> d={d} < {d_bound}"

res = wiener_attack(e, N)
print("CLASSICAL ",res)
#d_recovered, p_recovered, q_recovered = wiener_attack(e, N)

#if d_recovered:
 #   print(f"Found secret exponent d: {d_recovered}") # got 3 instead of 5 i wonder y
 #   print(f"Prime factors p: {p_recovered}, q: {q_recovered}")
    #assert d == d_recovered, "not as original"

#else:
 #   print("d not found")

#bits = 50

#p = 107288913396737 #random_prime(2**bits)
#q = 448210034652803 #random_prime(2**bits)
#N = p * q
#phi_N = (p - 1) * (q - 1)
#d = 3
#e = inverse_mod(d, phi_N)
#assert (e * d) % phi_N == 1, f"e * d does not satisfy ed ≡ 1 (mod φ(N))  -> {(e * d) % phi_N } "
#assert d < int((1/3) * N**(1/4)), f"d does not satisfy the condition d < 1/3 * N^(1/4) d={d}, 1/3 * N^(1/4) = { int((1/3) * N**(1/4))} "


# Output the values
#print(f"p: {p}")
#print(f"q: {q}")
#print(f"N: {N}")
#print(f"φ(N): {phi_N}")
#print(f"e: {e}")
#print(f"d: {d}")

#resL =  wiener_attack_lattice(e, N)
#print("LATTICE",resL)