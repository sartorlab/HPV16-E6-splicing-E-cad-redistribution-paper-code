#!/usr/bin/env python3
"""Figure 7E: Overall and recurrence-free survival in HPV16-positive patients.

Input tables: figure7/data/survival_n230.csv
Scientific code: figure7/scripts/analyze_survival_and_groups.R; figure7/scripts/render_survival_and_groups.py
Method: HPV16 population-mean score split; continuous-score Cox adjustment for age, stage and cohort.
Output: Figure7E.pdf, plus its underlying statistical tables.

Example: python panels/Figure7E.py --output build/Figure7E
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="7E")
