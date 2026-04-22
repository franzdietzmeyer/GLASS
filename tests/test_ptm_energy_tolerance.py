import sys
from pathlib import Path

import pandas as pd

# Add analysis/ so we can import ptm_analysis
ANALYSIS = Path(__file__).resolve().parent.parent / "analysis"
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))

from ptm_analysis import PTMAnalyzer  # noqa: E402


def test_ptm_energy_filter_tolerance_keeps_within_pct():
    a = PTMAnalyzer(debug=False)
    df = pd.DataFrame(
        {
            "description": ["a_10_0001", "a_10_0002", "a_10_0003"],
            "native_total_energy": [-100.0, -100.0, -100.0],
            "final_total_energy": [-95.0, -89.0, -100.0],  # +5%, +11%, 0%
            "position": [10, 10, 10],
        }
    )

    # Allow up to 10% worse: keep delta<=10 (abs(pre)*0.1) => first & third remain.
    out = a._filter_ptm_rows_worse_post_energy(df, allow_worse_scores_pct=10.0)
    assert len(out) == 2
    assert set(out["description"]) == {"a_10_0001", "a_10_0003"}


def test_ptm_energy_filter_strict_excludes_any_worse():
    a = PTMAnalyzer(debug=False)
    df = pd.DataFrame(
        {
            "description": ["b_20_0001", "b_20_0002"],
            "native_total_energy": [-100.0, -100.0],
            "final_total_energy": [-99.0, -101.0],  # worse, better
            "position": [20, 20],
        }
    )
    out = a._filter_ptm_rows_worse_post_energy(df, allow_worse_scores_pct=0.0)
    assert len(out) == 1
    assert out.iloc[0]["description"] == "b_20_0002"

