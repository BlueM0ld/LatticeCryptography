import os
import uuid
import shutil
import subprocess


def update_filename(dir):
    """"
        This will update the filename between runs for lattice instances. Keeping run executions separate
    """

    dir_path = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), dir)
    print(f"Directory to be checked: {dir_path}")

    if not os.path.exists(dir_path):
        print(f"Directory '{dir_path}' does not exist.")
        return

    for root, dirs, files in os.walk(dir_path):
        for filename in files:
            if (dir == "../result"):
                # Check if the file has the .npz extension
                res_dir(root, filename)
            if (dir == "../src"):
                # Check if the file has the .sage.py extension
                src_dir(root, filename)


def src_dir(root, filename):
    """"
    Pre execution script need to generate py files for execution likely overkill with rechecking but will do ...
    """
    if filename.endswith('.sage'):
        sage_filepath = os.path.join(root, filename)
        sage_py_filepath = sage_filepath + '.py'

        # Execute the .sage file to create .sage.py file
        subprocess.run(['sage', '--preparse', sage_filepath], check=True)

        py_filepath = sage_filepath.replace('.sage', '') + '.py'
        os.rename(sage_py_filepath, py_filepath)


def res_dir(root, filename):
    """
    Adds uuid to the filename, helper function to update_file_name
    """
    if filename == "reduced_gso_norms.npz":
        ran_str = str(uuid.uuid4())

        new_filename = f"reduced_gso_norms_{ran_str}.npz"

        old_filepath = os.path.join(root, filename)
        new_filepath = os.path.join(root, new_filename)

        os.rename(old_filepath, new_filepath)


def delete_folder():
    """
    Deletes the result folder ideally should execute at a sub directory level to give more control
    """
    res_dir = os.path.join(os.path.abspath(
        os.path.dirname(__file__)), "../result")

    if os.path.exists(res_dir):
        shutil.rmtree(res_dir)
        print(f"Deleted folder: {res_dir}")
    else:
        print(f"Folder does not exist: {res_dir}, or files are read-only?")
