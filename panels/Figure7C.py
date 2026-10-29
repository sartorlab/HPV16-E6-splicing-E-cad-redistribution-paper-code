#!/usr/bin/env python3
"""Figure 7C: E6 score comparison between IMU and KRT subtypes.

Input tables: figure7/data/subtype_integration_n230.csv
Scientific code: figure7/scripts/analyze_survival_and_groups.R; figure7/scripts/render_survival_and_groups.py
Method: Two-sided Student t-test on the annotated patient subset.
Output: Figure7C.pdf, plus its underlying statistical tables.

Example: python panels/Figure7C.py --output build/Figure7C
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="7C")
