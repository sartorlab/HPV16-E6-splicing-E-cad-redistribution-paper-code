#!/usr/bin/env python3
"""Figure S7C: Volcano plot: E6*I versus PCD.

Input tables: cell_line/data/de/res_03_E6IvsPCD.tsv
Scientific code: cell_line/scripts/04_volcano_plots.R
Method: Supplied DESeq2 results; adjusted p < 0.05 and absolute log2 fold change > 1.
Output: FigureS7C.pdf, plus its underlying statistical tables.

Example: python panels/FigureS7C.py --output build/FigureS7C
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="S7C")
