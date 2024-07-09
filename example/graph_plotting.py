from sage.all import *
import matplotlib.pyplot as plt
from fpylll import IntegerMatrix, GSO


def plot_gso(original_log_gso_norms, reduced_log_gso_norms, output_file="gso_plot.png"):
    plt.figure(figsize=(14, 6))

    # subplot 1 - M Orginal GSO norms
    plt.subplot(1, 2, 1)
    for i, vec in enumerate(original_log_gso_norms):
        plt.plot(range(len(vec)), vec, label=f"Original Vector {i+1}")
    plt.ylabel("log Gram-Schmidt Norms")
    plt.title("Original GSO Norms")
    plt.legend()

    # subplot 2 - Reduced GSO norms
    plt.subplot(1, 2, 2)
    for i, vec in enumerate(reduced_log_gso_norms):
        plt.plot(range(len(vec)), vec,
                 label=f"Reduced Vector {i+1}", linestyle='--')
    plt.ylabel("log Gram-Schmidt Norms")
    plt.title("Reduced GSO Norms - LLL Reduction")
    plt.legend()

    plt.savefig(output_file)
    plt.show()
    plt.close()


def compute_and_plot_gso(M, output_file="gso_plot"):
    fpylll_matrix = convert_to_fpylll(M)
    print("fpylll matrix", fpylll_matrix)
    M_gso = GSO.Mat(fpylll_matrix)
    print("GSO norms", M_gso)
    M_gso.update_gso()
    square_gso_norms = M_gso.r()
    original_log_gso_norms = [(log(square_gso_norm, 2)/2)
                              for square_gso_norm in square_gso_norms]

    original_log_gso_norms = [N(norm) for norm in original_log_gso_norms]

    print("Performing LLL reduction...")
    reducedM = M.LLL()

    reducedM_fpylll_matrix = convert_to_fpylll(reducedM)
    print("fpylll matrix", reducedM_fpylll_matrix)
    reducedM_gso = GSO.Mat(reducedM_fpylll_matrix)
    print("GSO norms", reducedM_gso)
    reducedM_gso.update_gso()
    red_square_gso_norms = reducedM_gso.r()
    reduced_log_gso_norms = [RR(log(red_square_gso_norm, 2)/2)
                             for red_square_gso_norm in red_square_gso_norms]
    # reduced_log_gso_norms = [N(norm) for norm in reduced_log_gso_norms]

    plot_gso([original_log_gso_norms], [reduced_log_gso_norms],
             "res/" + output_file + ".png")

    print("Reduced Matrix generated:")
    print(reducedM.str(rep_mapping=lambda x: str(x.n(digits=3))))
    print("COL:", reducedM.ncols())
    print("ROWS:", reducedM.nrows())
    print("DIMENSIONS:", reducedM.dimensions())

    return reducedM


def convert_to_fpylll(mat):
    # There's an issue with infinite values - this is arbitary as can't find a limit!
    # TODO: still running into issue need to think what to do here
    # max_val = 10**20
    # min_val = -10**20
    # safe_mat = mat.apply_map(lambda x: min(max_val, max(min_val, x)) if x not in [
    #                        float('inf'), float('-inf')] else max_val if x == float('inf') else min_val)
    return IntegerMatrix.from_matrix(mat)
