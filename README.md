# Lattice Cryptanalysis Project

This project is an artifact for a master's level thesis. It explores cryptanalytic attacks on a relaxed RSA cryptosystem using lattice-based methods, which can be executed from the terminal.

![Terminal screenshot](misc/image-2.png)

## Table of Contents

1. [Getting Started](#getting-started)
   - [Prerequisites](#prerequisites)
2. [Usage](#usage)
3. [License](#license)
4. [Acknowledgments](#acknowledgments)
5. [References](#references)


## Getting Started

These are the steps to get this project running on your machine. Installation steps may vary slightly based on your operating system (Windows, Mac, etc.).

### Prerequisites

To run this project, you will need the following:

-  **Python**: This project uses Python 3.10.12

-  **SageMath**: Install SageMath via Conda, you can install by other means but I'm not sure if affect the setup. Follow instructions here under the Not for development
[SageMath Conda](https://doc.sagemath.org/html/en/installation/conda.html#installing-all-of-sagemath-from-conda-not-for-development)

### Installation

1. **Clone the repository**:
   ```
   git clone git@github.com:BlueM0ld/LatticeCryptography.git
   cd lattice-cryptanalysis-project
   ```

2. **Install requirements**:
   ```
   pip install -r requirements.txt
   ```


### Running

To run the script to execute attacks following these steps, these are executed at the root of the project. This has a pre-execution script that will ensure setup if you want to run the experiments.

   ```
    cd src
    sage run.sage
   ```

To run the set experiments you will need to uncomment the tests and run the experiments one at a time. After each run delete result folder. Or you can move the file created to another space. This is because data will be extended. Replace {experiment_file_name} will the file that you want to execute.


   ```
    cd experiments
    sage {experiment_file_name}.sage
   ```

 
## License

Distributed under the MIT License. See `LICENSE` for more information.

## References

Here are the papers that we're used in the implementation of these attacks.

