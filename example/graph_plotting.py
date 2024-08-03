from sage.all import *
import matplotlib.pyplot as plt
from fpylll import IntegerMatrix, GSO

debug_flag = False


# Plots original/reduced GSO norms
def plot_gso(original_log_gso_norms=None, reduced_log_gso_norms=None, output_file="gso_plot.png"):

    if original_log_gso_norms is None and reduced_log_gso_norms is None:
        print("No valid GSO norms to plot.")
        return

    plt.figure(figsize=(14, 6))

    # subplot 1 - M Original GSO norms
    if original_log_gso_norms is not None:
        plt.subplot(1, 2, 1)
        for i, vec in enumerate(original_log_gso_norms):
            plt.plot(range(len(vec)), vec, label=f"Original Vector {i+1}")
        plt.ylabel("log Gram-Schmidt Norms")
        plt.grid(True)
        plt.title("Original GSO Norms")
        plt.legend()

    # subplot 2 - Reduced GSO norms
    if reduced_log_gso_norms is not None:
        plt.subplot(1, 2, 2)
        for i, vec in enumerate(reduced_log_gso_norms):
            plt.plot(range(len(vec)), vec,
                     label=f"Reduced Vector {i+1}", linestyle='--')
        plt.ylabel("log Gram-Schmidt Norms")
        plt.grid(True)
        plt.title("Reduced GSO Norms - LLL Reduction")
        plt.legend()

    plt.savefig(output_file)
    plt.show()
    plt.close()


# run LLL, get gso norms and call plot
def compute_and_plot_gso(M, output_file="gso_plot", debug=debug_flag):

    if debug:
        print("Matrix generated:")
        print(M.str(rep_mapping=lambda x: str(x.n(digits=3))))

    original_log_gso_norms = None
    reduced_log_gso_norms = None

    try:
        original_log_gso_norms = compute_log_gso_norms(M)
    except Exception as e:
        print(f"Error computing original GSO norms: {e}")

    reducedM = M.LLL()
    try:
        reduced_log_gso_norms = compute_log_gso_norms(reducedM)
    except Exception as e:
        print(f"Error computing reduced GSO norms: {e}")

    if original_log_gso_norms is not None or reduced_log_gso_norms is not None:
        plot_gso(
            original_log_gso_norms=[
                original_log_gso_norms] if original_log_gso_norms is not None else None,
            reduced_log_gso_norms=[
                reduced_log_gso_norms] if reduced_log_gso_norms is not None else None,
            output_file="res/" + output_file + ".png"
        )

    if debug and reducedM is not None:
        print("Reduced Matrix generated:")
        print(reducedM.str(rep_mapping=lambda x: str(x.n(digits=3))))
        print("COL:", reducedM.ncols())
        print("ROWS:", reducedM.nrows())
        print("DIMENSIONS:", reducedM.dimensions())

    return reducedM


# computes the log GSO norms for a given matrix
def compute_log_gso_norms(M, debug=debug_flag):

    fpylll_matrix = convert_to_fpylll(M)
    M_gso = GSO.Mat(fpylll_matrix)
    M_gso.update_gso()
    square_gso_norms = M_gso.r()
    if debug:
        print("square gso norm", square_gso_norms)
    log_gso_norms = [RR((log(norm, 2) / 2) if norm !=
                        0 else float('0')) for norm in square_gso_norms]
    return log_gso_norms


# converts a sage matrix to a fpylll matrix
def convert_to_fpylll(mat):
    return IntegerMatrix.from_matrix(mat)
