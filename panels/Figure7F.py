#!/usr/bin/env python3
"""Figure 7F: The selected 16-association network in Bailey's figure style.

Input tables: figure7/data/network_n230.csv; figure7/data/network_candidates.csv
Scientific code: figure7/scripts/analyze_network.R; figure7/scripts/render_network.py
Method: Eighteen fixed candidate associations adjusted for sequencing batch; nominal p < 0.05.
Output: Figure7F.pdf, plus its underlying statistical tables.

Example: python panels/Figure7F.py --output build/Figure7F
Common analysis and plotting functions are shared across related panels.
"""

from pathlib import Path
import sys

# Locate shared code independently of the caller's current working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from reproduce import main

if __name__ == "__main__":
    main(panel="7F")
