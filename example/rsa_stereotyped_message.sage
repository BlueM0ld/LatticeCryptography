from sage.all import *

# Define the RSA encryption parameters
N = random_prime(2**150) * random_prime(2**150)  # RSA modulus N
message = Integer('thepasswordfortodayisswordfish', base=35)  # Message to encrypt

# Encrypt the message using RSA encryption (with e=3)
c = message^3 % N

# Define the known and unknown parts of the message
known_prefix = 'thepasswordfortodayis000000000'
unknown_part = 'xxxxxxxxx'

# Convert known and unknown parts to integers in base 35
a = Integer(known_prefix, base=35)
X = Integer(unknown_part, base=35)

# Construct the matrix M for LLL reduction
M = matrix([
    [X^3, 3*X^2*a, 3*X*a^2, a^3 - c],
    [0, N*X^2, 0, 0],
    [0, 0, N*X, 0],
    [0, 0, 0, N]
])

# Perform LLL reduction on M to obtain the reduced basis B
B = M.LLL()

# Extract the first polynomial Q from the reduced basis B
Q = B[0][0]*x^3/X^3 + B[0][1]*x^2/X^2 + B[0][2]*x/X + B[0][3]

# Find the roots of Q(x) over the ring of integers (ZZ)
recovered_root = Q.roots(ring=ZZ)[0][0]

# Convert the recovered root to a base-35 string to recover the message
recovered_message = Integer(recovered_root).str(base=35)

# Print the recovered message
#print("B:", B)
print("Q:", Q)
print("Recovered root:", recovered_root)
print("Recovered Message:", recovered_message)
