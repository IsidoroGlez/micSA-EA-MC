# PROGRAM: Microcanonical Simulated Annealing Monte Carlo Simulation on GPU

This folder contains the source code of the simulation program for studing the 3 dimensional Edwards-Andersiong spin glass model of liner size L. The simulation use the microcanonical Simulated Annealing algorithm presented in [DOI]. It is written in **C/CUDA** and is designed to run on **NVIDIA GPUs** with CUDA support.

---

## Folder Structure

├── Makefile		# Build system
├── include/ 		# Header files
├── input/ 		# Input files (beta, configuration, LUT, time lists)
├── src/ 		# Source code (C and CUDA)
├── obj/ 		# Object files (created after compilation)
└── bin/ 		# Executable binary (created after compilation)

---

## Compilation

To build the simulation binary, simply run:

```bash

make

```

This will generate a file named like:

```bash

bin/EA_CUBE_L160_NB1_NBPRE10

```

The values L (system linear size) , NBETAS (number of temperatures ---the actual version in develop for use NBETAS=1), and NBITSPREBUSQUEDAS (number of bits use for searching in the Look-Up tables for generating random numbers) are defined at the top of the Makefile and passed as macros during compilation.

Make sure your system has:

 - nvcc from the NVIDIA CUDA Toolkit

 - A standard g++ compiler

 - The Makefile is preconfigured to support CUDA architecture sm_90. If needed, modify the NVCC_FLAGS for other GPU generations.

---

## Running the Program

The compiled binary requires between 7 and 9 arguments (some optional). To see usage details:

```bash

./bin/EA_CUBE_L160_NB1_NBPRE10 -h

```

### Example call:

```bash

./bin/EA_CUBE_L160_NB1_NBPRE10 \
    0 10 input/BETAS/beta.in input/input.in \
    input/LUT/LUT_for_PRNG_nbits10_NB1.bin input/measures_times.in \
    0 0 input/sample_list.txt
```

### Argument list:

| Pos | Argument       | Type       | Description                                         |
| --- | -------------- | ---------- | --------------------------------------------------- |
| 1   | `isample`      | `int`      | Sample index                                        |
| 2   | `nbits`        | `int`      | Number of Subsamples to run in parallel             |
| 3   | `beta.dat`     | `string`   | Path to `beta.in` file                              |
| 4   | `input.in`     | `string`   | Path to the main simulation input file              |
| 5   | `LUT`          | `string`   | Path to precomputed lookup table (`.bin`)           |
| 6   | `list_times`   | `string`   | File with list of measure times                     |
| 7   | `device`       | `int`      | CUDA device number                                  |
| 8   | `max_time`     | `unsigned` | *(Optional)* Max simulation time (0 disables limit) |
| 9   | `list_samples` | `string`   | *(Optional)* File with list of sub-samples          |


### Preparing Input Files

#### Lookup Table (LUT)

The simulation requires a precomputed Lookup Table (LUT) to speed up the random number generation. This table is provided as a binary .bin file and must be generated before running the main program.

To build the LUT, you must first compile the generation code using the provided script `input/LUT/compiler.sh`.

This script takes two arguments:

```bash

./compiler.sh NBITSPREBUSQUEDA NBETAS

```

- NBETAS: Number of distinct inverse temperature values (must match the count in beta.in ---this version only support NBETAS=1).

- NBITSPREBUSQUEDA: Number of bits used for pre-search (typically 10 or 11).

Once compiled, the LUT must be generated using the executable creatin_nbitsLUT{NBITSPREBUSQUEDA}_NB{NBETAS} that is created in the process. This executable expects a single argument:

```bash

./creatin_nbitsLUT10_NB01 path/to/beta.in

```

This will produce the final .bin file to be passed as the LUT argument to the main simulation binary.

#### Temperature file

The file `input/BETAS/betas.in` contain the value of the inverse temperature used in the simulation.

The actual version is develop and tested only for use a single temperature, then NBETAS=1.

---

## Notes

You can modify simulation parameters such as L, NBETAS (not recommended in the actual version), or NUMBITSPREBUSQUEDAS at the top of the Makefile.

The simulation uses precomputed lookup tables found under input/LUT/.

All timing, energy measurement, and sampling parameters are specified in input/*.in files.

Compilation will create the obj/ and bin/ directories automatically if needed.

---

## Cleaning Up

To remove object files and the compiled binary, run:

```bash

make clean

```

---

## Requirements

- CUDA Toolkit (tested with nvcc)

- C++ compiler (g++)

- GPU with appropriate compute capability (e.g., sm_90)

- A Linux or WSL environment is strongly recommended. Windows native support is not tested.

---

## License

This program is released under the MIT License. See ../LICENSE for details.
