# Threshold simulation of Hyperbolic codes

#Careful: the code has still some bugs and still needs to be tested thoroughly!

This repository generates hyperbolic codes and uses the autotune library https://github.com/adamcw/autotune) to perform a threshold simulation. Errors can happen at any time, such that syndrome extraction is not assumed to be ideal. Even if the reader is not interested in hyperbolic codes, this repository is might still be relevant because we provide a wrapper around the autotune library such that only lists of syndromes need to be provided for the simulation.



## Generation of Hyperbolic code patches
The files 
The code in this repository can be used to generate patches of hyperbolic codes. One can generate arbitrary sizes and types of lattices definde by the Schäffli-symbols {r,s}. The code generates the patches layer by layer, starting at the center and growing outwards. Then rough and smooth edges are generated on the boundary such that logical qubits can be encoded.

## Wrapper for autotune


# Installation and usage
## Building the library
To install this library with tests cmake version 3.11 is required. Furthermore, we require that the autotune is properly installed on the system. Inside the `CMakeLists.txt` file the folder location of the autotune library can be edited.

It is advised to build out of source. Thus, you should create a build directory:

```
mkdir build
cd build
```

Now the current folder is the newly created build folder we can instruct cmake to generate the build system:

```
cmake ../
```

Finally we can everything using Gnu make:

```
make all
```

## Running the program

The program can be run from the build directory using the command:
```
Usage:
        ./src/main      [--OPTIONS]


Options:
        --layers        number of layers and thus the size of the code
        --k             number of neighbors each surface has
        --r             number of neighbors each vertex has
        --t_check       The number of timesteps between each round of perfect stabilizer measurements
        --big_t_max     The maximum value of big_t
        --probability   The probability of a random error
        --boot          Flag if boot up phase should be performed to optimize t_check
        --seed0         The first random seed
        --seed1         The second random seed
        --help          Prints this help message

Standard parameters:
        ./src/main --layers 1 -k 4 -r 4 --t_check 10 --big_t_max 20000 --probability 0.001 --seed0 42 --seed1 42
```