from sage.all import *
import matplotlib.pyplot as plt
from fpylll import IntegerMatrix, GSO
import os
import numpy as np
import time
from fpylll.tools.quality import basis_quality
import sys


def get_basis_quality(M, debug=False):
    bq = basis_quality(M)
    listBq = [(k, v) for k, v in bq.items()]
    if debug:
        print("basis quality: ", listBq)
    return bq


def compute_and_plot_gso(M, folder, reduction):

    # print(reduction)

    print(f"Matrix Dimension: {M.dimensions()} ")

    # Check for result folder

    # res_dir = os.path.join("result", folder, reduction)
    res_dir = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), "../result", folder, reduction)
    os.makedirs(res_dir, exist_ok=True)

    reduced_log_gso_norms = None

    start_time, reducedM, end_time = reduction_method(M, reduction)

    elapsed_time = end_time - start_time
    print("elapsed time: ", elapsed_time)

    try:
        compute_log_gso_norms(reducedM, res_dir)
    except Exception as e:
        print(f"Error computing reduced GSO norms: {e}")

    return reducedM


def scale_down(M, gcd):
    if (gcd > 1):
        T = M.apply_map(lambda x: x // gcd)
        T.change_ring(ZZ)
        return T
    return M


def reduction_method(M, reduction):
    if reduction == "LLL":
        print("Performing LLL reduction...")
        start_time = time.time()
        reducedM = M.LLL()  # algorithm="fpLLL:heuristic", fp="rr", prec=5050)
        end_time = time.time()
    if reduction == "BKZ":
        print("Performing BKZ reduction...")
        start_time = time.time()
        reducedM = M.BKZ(precision=8000)
        end_time = time.time()
    if reduction == "BKZ40":
        print("Performing BKZ 40 reduction...")
        start_time = time.time()
        reducedM = M.BKZ(block_size=40, fp="rr",
                         precision=8000)
        end_time = time.time()
    if reduction == "BKZ60":
        print("Performing BKZ 60 reduction...")
        start_time = time.time()
        reducedM = M.BKZ(block_size=60, fp="rr",
                         precision=8000)
        end_time = time.time()
    return start_time, reducedM, end_time


# returns the log GSO norms for a given matrix
def compute_log_gso_norms(M, res_dir):

    fpl_M = convert_to_fpylll(M)

    M_gso = GSO.Mat(fpl_M)
    M_gso.update_gso()

    square_gso_norms = M_gso.r()
    # print("square norms", square_gso_norms)
    reduced_log_gso_norms = []
    bas_quality = []
    try:
        reduced_log_gso_norms = [(RR((log(norm, 2))) if norm !=
                                  0 else float('0')) for norm in square_gso_norms]
        bas_quality = get_basis_quality(M_gso)

        save_data([[M_gso.d], [reduced_log_gso_norms], bas_quality],
                  res_dir, "reduced_gso_norms.npz")

    except Exception as e:
        print(
            f"Error computing data for the GSO norms, exiting to save integrity of previous results: {e}")
        # Exit Gracefully!
        return

    return  # reduced_log_gso_norms, M_gso.d, bas_quality


def convert_to_fpylll(mat):
    return IntegerMatrix.from_matrix(mat)


def save_data(data, dir, filename):
    path_file = os.path.join(dir, filename)

    np_data = np.array(data, dtype=object)

    if os.path.exists(path_file):
        existing_data = np.load(path_file, allow_pickle=True)
        existing_np_data = existing_data['arr_0']

        # Handle varying lengths and concatenate
        # Convert to lists to handle varying lengths properly
        combined_data = list(existing_np_data) + list(np_data)
        np.savez(path_file, arr_0=np.array(combined_data, dtype=object))

    else:
        np.savez(path_file, np_data,  allow_pickle=True)


def calculate_gso(red_norms):
    # print("red_norms", red_norms)
    pred_gsa = []

    for norms in red_norms:
        b1 = norms[0]  # First element of the sublist
        # Normalize each element
        normalized_norms = [norm-b1 for norm in norms]
        pred_gsa.append(normalized_norms)

    print(pred_gsa)


def generate_graphs(attack, reduction, filename):
    path_to_file = folder_check(attack, reduction, filename)

    dimensions, red_norms, rhf = extract_data(path_to_file)
    print("rhf", rhf)
    # pred_gso = calculate_gso(red_norms)
    plot_graph(attack, reduction, dimensions, red_norms)


def plot_graph(attack, reduction, dimensions, red_norms):
    for norms, dim in zip(red_norms, dimensions):
        plt.figure(figsize=(7, 6))
        plt.plot(norms, label=f'Dimension {dim}')
        plt.xlabel('Index')
        plt.ylabel('Log GSO Norm')
        plt.title(f'Log GSO Norms after {reduction} Reduction')
        plt.grid()
        plt.legend()
        plt.show()
        plt.close()

#        filename = f"{attack}_{reduction}_{dim}.png"
#        output_file = os.path.join(res_dir, filename)
#       plt.savefig(output_file)


def combine_gso_norms(attack, reduction, filename):
    path_to_file = folder_check(attack, reduction, filename)

    dimensions, red_norms, rhf = extract_data(path_to_file)
    iter = list(range(0, len(dimensions)))

    # plot_red_gso_norms_iter(iter, red_norms, reduction)
    # plot_bas_q_iter(iter, rhf)
    plot_iter_combined(iter, red_norms, reduction, rhf)


def extract_data(path_to_file):
    reduced_log_gso_norms = np.load(path_to_file, allow_pickle=True)

    data = reduced_log_gso_norms['arr_0']

    # X values corresponding to each dataset
    dimensions = [data[i-1][0] for i in range(1, len(data), 3)]
    red_norms = []

    for i in range(1, len(data), 3):
        y = data[i][0]
        # Handle infinities: Replace np.inf with a large number for plotting
        y = [val if val != np.inf else max(y) + 1 for val in y]
        red_norms.append(y)

    rhf_list = []

    for i in range(2, len(data), 3):
        rhf = data[i]['rhf']
        rhf_list.append(rhf)
    return dimensions, red_norms, rhf_list


def folder_check(attack_type, reduction_type, filename):
    result_directory = os.path.join("../result", attack_type, reduction_type)
    os.makedirs(result_directory, exist_ok=True)

    if filename is None:
        filename = "reduced_gso_norms"

    file_path = os.path.join(result_directory, f"{filename}.npz")

    # Check if the file exists
    if not os.path.isfile(file_path):
        print(
            f"File does not exist: {file_path}, please check or run a reduction first!")
        sys.exit(1)  # Exit Gracefully!

    # Return the path to the file if it exists
    return file_path


def plot_red_gso_norms_iter(iter, red_norms, reduction):
    plt.figure(figsize=(10, 6))
    for norms, i in zip(red_norms, iter):
        plt.plot(norms, label=f'Interation {i}')
    plt.ylabel('Log GSO Norm')
    plt.title(f'Log GSO Norms after {reduction} Reduction')
   # plt.legend()
    plt.grid()
    plt.show()
    plt.close()


def plot_iter_combined(iter, red_norms, reduction, rhf):
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))

    # First subplot: Log GSO Norms
    for norms, i in zip(red_norms, iter):
        axs[0].plot(norms, label=f'Iteration {i}')
    axs[0].set_ylabel('Log GSO Norm')
    axs[0].set_title(f'Log GSO Norms after {reduction} Reduction')
    axs[0].grid(True)
    # axs[0].legend() # Uncomment to show legend

    # Second subplot: Root Hermite Factor
    axs[1].plot(iter, rhf, marker='o', label='Root Hermite Factor')
    axs[1].set_xlabel('Iteration')
    axs[1].set_ylabel('Root Hermite Factor')
    axs[1].set_title('Root Hermite Factor per Iteration')
    axs[1].grid(True)
    # axs[1].legend() # Uncomment to show legend

    plt.tight_layout()
    plt.show()
    plt.close()


def plot_bas_q_iter(iter, rhf):

    plt.figure(figsize=(10, 6))
    plt.plot(iter, rhf, marker='o', label=f'Iteration {iter}')
    plt.xlabel('Iteration')
    plt.ylabel('Root Hermite Factor')
    plt.title('Root Hermite Factor per Iteration')
    # plt.legend()
    plt.grid()
    plt.show()
    plt.close()


def plot_red_gso_norms_dimension(dimensions, red_norms, reduction):
    plt.figure(figsize=(10, 6))
    for norms, dim in zip(red_norms, dimensions):
        plt.plot(norms, label=f'Dimension {dim}')
    plt.ylabel('Log GSO Norm')
    plt.title(f'Log GSO Norms after {reduction} Reduction')
   # plt.legend()
    plt.grid()
    plt.show()
    plt.close()


def plot_bas_q_dimension(dimensions, rhf):
    plt.figure(figsize=(10, 6))
    plt.plot(dimensions, rhf, marker='o', label=f'Dimension {dimensions}')
    plt.xlabel('Dimension')
    plt.ylabel('Root Hermite Factor')
    plt.title('Root Hermite Factor vs. Dimension')
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()
