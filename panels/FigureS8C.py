#!/usr/bin/env python3
"""Figure S8C: GO-BP, Hallmark and KEGG enrichment: E6*I versus PCD.

Input tables: cell_line/data/ora/all_ORA_results.tsv
Scientific code: cell_line/scripts/05_ora_plots.R
Method: Supplied ORA results; top 10 significant terms per direction and database.
Output: FigureS8C.pdf, plus its underlying statistical tables.

Example: python panels/FigureS8C.py --output build/FigureS8C
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="S8C")
