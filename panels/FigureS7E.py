#!/usr/bin/env python3
"""Figure S7E: Volcano plot: FLm versus PCD.

Input tables: cell_line/data/de/res_06_FLmvsPCD.tsv
Scientific code: cell_line/scripts/04_volcano_plots.R
Method: Supplied DESeq2 results; adjusted p < 0.05 and absolute log2 fold change > 1.
Output: FigureS7E.pdf, plus its underlying statistical tables.

Example: python panels/FigureS7E.py --output build/FigureS7E
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="S7E")
