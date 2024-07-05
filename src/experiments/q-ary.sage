import sys
import os

src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(src_path)

import sage.all as sage
from fpylll import IntegerMatrix, LLL, GSO
from graph_plotting import compute_and_plot_gso

# Set random seed for reproducibility
sage.set_random_seed(1337)

def generate_qary_lattice(size):
    A = IntegerMatrix.random(size, "qary", k=size//2, bits=40)
    A_sage = sage.matrix(A)
    return A_sage

print("Generating q-ary lattices...")

# defining the range of lattice sizes
# 50 -100 results came out as expected not much variation okay
# 100-150 results came out as expected not much variation again
# sanity check done
start_size = 100
end_size = 150

for size in range(start_size, end_size + 1):
    lattice = generate_qary_lattice(size)
    print(f"Lattice generated with size {size}:")
    print(lattice)
    compute_and_plot_gso(lattice, "q-ary",f"q-ary_{size}")

print("All q-ary lattices generated.")
