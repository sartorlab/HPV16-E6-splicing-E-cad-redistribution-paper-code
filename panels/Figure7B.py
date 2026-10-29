#!/usr/bin/env python3
"""Figure 7B: Overall and recurrence-free survival in all HPV-positive patients.

Input tables: figure7/data/survival_n230.csv
Scientific code: figure7/scripts/analyze_survival_and_groups.R; figure7/scripts/render_survival_and_groups.py
Method: Population-mean score split; continuous-score Cox adjustment for age, stage and cohort.
Output: Figure7B.pdf, plus its underlying statistical tables.

Example: python panels/Figure7B.py --output build/Figure7B
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="7B")
