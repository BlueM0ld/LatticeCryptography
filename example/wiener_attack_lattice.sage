from sage.all import *
from graph_plotting import compute_and_plot_gso, convert_to_fpylll 
from fpylll import IntegerMatrix, GSO


def wiener_attack_lattice(e, N, trued):
    # Define the polynomial
    PR.<u, x> = PolynomialRing(ZZ)
    f = u + x * N

    # Set up bounds based on the description
    X = floor(N**(1/4))  # Bound for k
    Y = floor(3 * N**(1/2))  # Bound for p + q - 1

    L = matrix(ZZ, [
        [e, int(sqrt(N))], 
        [N, 0]
    ])


    # Apply LLL algorithm to reduce the lattice basis
    print(f"Lattice basis matrix :\n{L}")

    L = compute_and_plot_gso(L, "wieners")


    u0, x0 = L[0]
    print(f"u0: {u0}, x0: {x0}")

    if gcd(u0, x0) != 1:
        print("The gcd of u0 and x0 is not 1, the result may not be correct.")

    d = x0 // int(sqrt(N))
    d = abs(d)
        
    if d == trued:
        print("Recovered private exponent d is correct.")
    return d


bits = 25

#p = 107288913396737 
#q = 448210034652803 
p = random_prime(2**bits)
q = random_prime(2**bits)
N = p * q
phi_N = (p - 1) * (q - 1)
d = 3
e = inverse_mod(d, phi_N)
assert (e * d) % phi_N == 1, f"e * d does not satisfy ed = 1 (mod phi(N))  -> {(e * d) % phi_N } "
assert d < int((1/3) * N**(1/4)), f"d does not satisfy the condition d < 1/3 * N^(1/4) d={d}, 1/3 * N^(1/4) = { int((1/3) * N**(1/4))} "


# Output the values
print(f"p: {p}")
print(f"q: {q}")
print(f"N: {N}")
print(f"phi(N): {phi_N}")
print(f"e: {e}")
print(f"d: {d}")

result = wiener_attack_lattice(e, N,d)
print(result)