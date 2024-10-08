# Coevolutionary Simulation

This repository contains the implementation of coevolutionary algorithms in both C++ and Python based on the model from 
["Coevolution of Node and Link States in Networks"](https://doi.org/10.1103/PhysRevE.103.052302).
The model explores the interplay between node and link states within complex networks.

## Contents
- `CPP_version/`: C++ implementation.
- `Python_version/`: Python implementation.
- `plot.py`: Script for generating relevant plots from data produced by the C++ version.

## Requirements (for Python version & plot)
- Python 3.x 
- Required Python packages:
  - numpy
  - matplotlib
  - numba


## Usage
### Python Version
1. Clone the repository:
    ```bash
    git clone https://github.com/MohiMozaffari/Coevolutionary-Simulation.git
    ```
2. Navigate to the Python version:
    ```bash
    cd Python_version
    ```
3. Install the required packages:
    ```bash
    pip install numpy matplotlib numba
    ```
4. Run the Python version of the model:
    ```bash
    python main.py
    ```
### C++ Version
1. Clone the repository:
    ```bash
    git clone https://github.com/MohiMozaffari/Coevolutionary-Simulation.git
    ```
2. Navigate to the C++ version:
    ```bash
    cd CPP_version
    ```
3. Compile the C++ code:
    ```bash
    g++ -std=c++11 main.cpp -o coevolutionary_model
    ```
4. Run the compiled C++ program:
    ```bash
    ./coevolutionary_model
    ```
