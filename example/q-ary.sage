from fpylll import IntegerMatrix, LLL, GSO
import sage.all as sage
import matplotlib.pyplot as plt
import time

def generate_lattice_qary(dim_size):
    sage.set_random_seed(1337)

    print("Generating a q-ary lattice...")

    q = sage.random_prime(2**30, false, 2**29)
    A = IntegerMatrix.random(dim_size, "qary", k=75, bits=30)
    M = GSO.Mat(A, float_type="d")

    print("Computing original GSO norms...")
    M.update_gso()

    print("Performing LLL reduction...")
    start_time = time.time()

    L = LLL.Reduction(M)
    L()

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("LLL reduction completed")
    print("Computing GSO norms after LLL reduction...")

    square_gso_norms_r = [M.get_r(i, i) for i in range(M.d)]
    reduced_log_gso_norms = [sage.RR(sage.log(norm, 2)/2) for norm in square_gso_norms_r]

    print("GSO norms after LLL reduction computed")

    return (dim_size, elapsed_time), reduced_log_gso_norms

# Initialize lists to store results
results = []
all_norms = []
dimensions = []
times = []

# Generate lattice data for specified dimensions
for i in range(100, 155, 5):
    res, norms = generate_lattice_qary(i)
    results.append(res)
    all_norms.append(norms)
    dimensions.append(res[0])
    times.append(res[1])
    print("RES: ", res)

# Plotting GSO norms and elapsed time per dimension
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={'width_ratios': [1.5, 1]})

# Plot GSO norms
for norms, dim in zip(all_norms, dimensions):
    ax1.plot(norms, label=f'Dimension {dim}')
ax1.set_xlabel('Index')
ax1.set_ylabel('Log GSO Norm')
ax1.set_title('Log GSO Norms after LLL Reduction')
ax1.legend()

# Plot elapsed time per dimension
ax2.plot(dimensions, times, linestyle='-')
ax2.set_xlabel('Dimension Size')
ax2.grid()
ax2.set_ylabel('Elapsed Time (s)')
ax2.set_title('Elapsed Time per Dimension')


plt.savefig("q-ary_lattice_reduction.png")
plt.show()
plt.close()
