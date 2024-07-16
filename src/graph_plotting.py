from sage.all import *
import matplotlib.pyplot as plt
from fpylll import IntegerMatrix, GSO
import os
import uuid
import numpy as np


def plot_gso(original_log_gso_norms=None, reduced_log_gso_norms=None, output_file=None, reduction=""):

    if output_file is None:
        output_file = f"gso_plot_{uuid.uuid4()}.png"

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
        plt.title("Original GSO Norms")
        plt.legend()

    # subplot 2 - Reduced GSO norms
    if reduced_log_gso_norms is not None:
        plt.subplot(1, 2, 2)
        for i, vec in enumerate(reduced_log_gso_norms):
            plt.plot(range(len(vec)), vec,
                     label=f"Reduced Vector {i+1}", linestyle='--')
        plt.ylabel("log Gram-Schmidt Norms")
        plt.title(f"Reduced GSO Norms - {reduction} Reduction")
        plt.legend()

    plt.savefig(output_file)
    # plt.show()
    plt.close()


def compute_and_plot_gso(M, folder, file_name=None, reduction="LLL"):

    # Check for result folder
    res_dir = os.path.join("result", folder, reduction)
    os.makedirs(res_dir, exist_ok=True)

    # Generate unique file name using UUID
    # TODO: need to update it to matcht expermiments rather than uuid
    if file_name is None:
        file_name = uuid.uuid4()

    output_file = os.path.join(res_dir, f"gso_plot_{uuid.uuid4()}.png")

    print("Matrix generated:")
    print(M.str(rep_mapping=lambda x: str(x.n(digits=3))))

    # M = scale_matrix(M)

    original_log_gso_norms = None
    reduced_log_gso_norms = None

    try:
        original_log_gso_norms = compute_log_gso_norms(M)
    except Exception as e:
        print(f"Error computing original GSO norms: {e}")

    if reduction == "LLL":
        reducedM = M.LLL()
    if reduction == "BKZ":
        reducedM = M.BKZ()
    try:
        reduced_log_gso_norms = compute_log_gso_norms(reducedM)
        save_data([reduced_log_gso_norms],  os.path.join(
            res_dir, "reduced_gso_norms.npy"))
    except Exception as e:
        print(f"Error computing reduced GSO norms: {e}")

    if original_log_gso_norms is not None or reduced_log_gso_norms is not None:
        plot_gso(
            original_log_gso_norms=[
                original_log_gso_norms] if original_log_gso_norms is not None else None,
            reduced_log_gso_norms=[
                reduced_log_gso_norms] if reduced_log_gso_norms is not None else None,
            output_file=output_file, reduction=reduction
        )

    if reducedM is not None:
        print("Reduced Matrix generated:")
        print(M.str(rep_mapping=lambda x: str(x.n(digits=3))))
        print("COL:", reducedM.ncols())
        print("ROWS:", reducedM.nrows())
        print("DIMENSIONS:", reducedM.dimensions())

    return reducedM

# scales the matrix by its GCD and changes the ring to ZZ.


def scale_matrix(M):

    gcd_M = M.gcd()
    M = M * (1/gcd_M)
    M = M.change_ring(ZZ)
    return M


# Computes the log GSO norms for a given matrix.
def compute_log_gso_norms(M):

    fpylll_matrix = convert_to_fpylll(M)
    print(fpylll_matrix)
    M_gso = GSO.Mat(fpylll_matrix)
    M_gso.update_gso()
    square_gso_norms = M_gso.r()
    print("square gso norm", square_gso_norms)
    log_gso_norms = [RR((log(norm, 2) / 2) if norm !=
                        0 else float('0')) for norm in square_gso_norms]
    return log_gso_norms


def convert_to_fpylll(mat):
    return IntegerMatrix.from_matrix(mat)


def save_data(data, filename):
    if os.path.exists(filename):
        existing_data = np.load(filename)
        updated_data = np.append(existing_data, data, axis=0)
    else:
        updated_data = data

    np.save(filename, updated_data)


def combine_gso_norms(attack, reduction):
    res_dir = os.path.join("result", attack, reduction)
    os.makedirs(res_dir, exist_ok=True)
    # Load existing data
    path_to_file = os.path.join(res_dir, "reduced_gso_norms.npy")
    reduced_log_gso_norms = np.load(path_to_file)
    print(reduced_log_gso_norms)
    if reduced_log_gso_norms is not None:
        plt.figure(figsize=(10, 6))
        for i, vec in enumerate(reduced_log_gso_norms):
            plt.plot(range(len(vec)), vec,
                     label=f"Reduced Vector {i+1}", linestyle='--')

        plt.xlabel('Index')
        plt.ylabel('log Gram-Schmidt Norms')
        plt.title('Reduced GSO Norms Comparison')
        plt.legend()

        output_file = os.path.join(res_dir, "collated_reduced_gso_norms.png")

        plt.savefig(output_file)
        plt.show()
        plt.close()
    else:
        print(f"File {path_to_file} does not exist.")
