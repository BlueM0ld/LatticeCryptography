from sage.all import *
import matplotlib.pyplot as plt
from fpylll import IntegerMatrix, GSO


def plot_gso(log_gso_norms, output_file="gso_plot.png"):
    plt.figure(figsize=(10, 6))
    for i, vec in enumerate(log_gso_norms):
        plt.plot(range(len(vec)), vec, label=f"Vector {i+1}")
    plt.ylabel("log Gram-Schmidt Norms")
    plt.title("LLL Reduction")
    plt.legend()
    plt.savefig(output_file)
    plt.close()


def compute_and_plot_gso(M, output_file="gso_plot"):
    fpylll_matrix = convert_to_fpylll(M)
    print("fpylll matrix", fpylll_matrix)
    M_gso = GSO.Mat(fpylll_matrix)
    print("GSO norms", M_gso)
    M_gso.update_gso()
    square_gso_norms = M_gso.r()
    log_gso_norms = [RR(log(square_gso_norm, 2)/2)
                     for square_gso_norm in square_gso_norms]

    plot_gso([log_gso_norms], "res/"+output_file+".png")

    reducedM = M.LLL()
    print("Reduced Matrix generated:")
    print(M.str(rep_mapping=lambda x: str(x.n(digits=3))))
    print("COL:", M.ncols())
    print("ROWS:", M.nrows())
    print("DIMENSIONS:", M.dimensions())

    return reducedM


def convert_to_fpylll(mat):
    # There's an issue with infinite values - this is arbitary as can't find a limit!
    max_val = 10**20
    safe_mat = mat.apply_map(lambda x: min(max_val, x)
                             if x != float('inf') else max_val)
    return IntegerMatrix.from_matrix(safe_mat)
