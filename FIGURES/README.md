# Figures for the Microcanonical Simulated Annealing Monte Carlo Project

This directory contains the **figure generation scripts and data** used in the article [DOI]. Each subdirectory (`fig1`, `fig2`, ...) corresponds to a specific figure, and includes:

- The Gnuplot script (`figX.gpt`)
- The final figure in PDF format (`figX.pdf`)
- The input data used for plotting (`DATA/`)
- Color palette definitions and gnuplot terminal configuration (`my_color_palete.gpt`, `APS_term.gpt`)

---

## How to generate a figure

To regenerate a figure, follow these steps:

1. **Go into the corresponding directory**, e.g.:

```bash

cd fig2

```

2. **Run the Gnuplot script**

```bash

gnuplot fig2.gpt

```

This will generate a LaTeX/EPS-based .tex and .eps files and compile it using pdflatex to produce the final fig2.pdf file.

---

## Requirements

- Gnuplot with support for LaTeX terminal (usually default).

- TeX distribution (e.g. TeX Live, MiKTeX, or similar) with pdflatex in your system path.

---

## Cleaning up

At the end of each Gnuplot script, a cleanup command (rm) removes temporary files like .log, .aux, and .tex. This works in Linux/macOS, but not directly in Windows. For running on Windows systems, substitute the Linux command `rm` with `del` in files figX.gpt.

## Folder structure

Each figure folder follows the same layout:

figX/
├── figX.gpt		# Gnuplot script
├── figX.pdf       	# Final generated figure
├── APS_term.gpt   	# Gnuplot formatting helpers
├── my_color_palete.gpt
└── DATA/		# Raw data used for the figure

## Notes

All figure files are version-controlled, but data files are excluded from Git history via .gitignore unless explicitly needed.

Output .pdf files are kept for convenience but can be regenerated as explained above.
