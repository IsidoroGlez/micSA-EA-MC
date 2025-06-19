# Microcanonical Simulated Annealing Monte Carlo Simulations in Spin Glasses

This repository contains the code and analysis tools used to perform **microcanonical Simulated Annealing** Monte Carlo simulations of spin glasses using CUDA-C.

The simulations and resulting figures correspond to [doi].

---

## Repository structure

├── PROGRAM/		# Source code and headers for the simulation (see PROGRAM/README.md)
├── FIGURES/ 		# Plotting scripts, data, and final figures (see FIGURES/README.md)
├── .gitignore 		# Exclusion rules for temporary and compiled files
└── README.md 		# This file

---

## What this project does

- Implements the CUDA-based microcanonical simulated Annealing Monte Carlo simulation algorithm for spin glass systems.
- Outputs raw data which is then post-processed to generate publication-quality figures.

---

## How to use this repository

1. **Compile and run the simulation**:
   - Go to `PROGRAM/` and follow the instructions in `PROGRAM/README.md`.

2. **Generate figures**:
   - Use `gnuplot` scripts located in `FIGURES/` to reproduce the plots used in the article. Read `FIGURES/README.md` for more information.

---

## License and citation

This project is licensed under the **MIT License**. See the LICENSE file for full license text.

If you use this code in a publication or derived work, please cite the associated article and/or this repository.

---

## Contact

Developed by M. Bernaschi, L.A. Fernandez, I. González-Adalid Pemartín, E. Marinari, V. Martín-Mayor, G. Parisi, F. Ricci-Tersenghi, J.J. Ruiz-Lorenzo and D. Yllanes

For questions, reach out at [isiglezadalid@gmail.com] or via GitHub issues.
